
using PedModule
using DataFrames
using Distributions

function Gibbs(A,x,b,nIter;outFreq=100)
    n = size(x,1)
    xMean = zeros(n)
    for iter = 1:nIter
        if iter%outFreq==0
            println("at sample: ",iter)
        end
        for i=1:n
            cVar = 1.0/A[i,i]
            cMean   = cVar*(b[i] - A[:,i]'x)[1,1] + x[i]
            x[i]    = randn()*sqrt(cVar) + cMean 
        end
        xMean += (x - xMean)/iter
    end
    return xMean
end

function Gibbs(A,x,b)
    n = size(x,1)
    for i=1:n
        cVar = 1.0/A[i,i]
        cMean   = cVar*(b[i] - A[:,i]'x)[1,1] + x[i]
        x[i]    = randn()*sqrt(cVar) + cMean 
    end
end

function GaussSeidel(A,x,b;tol=0.000001)
    n = size(x,1)
    for i=1:n
        x[i] = ((b[i] - A[:,i]'x)/A[i,i])[1,1] + x[i]
    end
    diff = sum((A*x-b).^2)
    iter = 0
    while ((diff/n > tol) & (iter<1000))
        iter += 1
        for i=1:n
            x[i] = ((b[i] - A[:,i]'x)/A[i,i])[1,1] + x[i]
        end
        diff = sum((A*x-b).^2)
        #println(iter," ",diff/n)
    end
    return x
end

function Jacobi(A,x,b,p;tol=0.000001)
    D       = diag(A)
    res     = A*x
    resid   = b-res
    tempSol = resid./D
    diff    = sum(resid.^2)
    n    = size(A,1)
    iter = 0
    while ((diff/n > tol) & (iter<1000))
        iter += 1
        x = p*tempSol + (1-p)*x
        res     = A*x
        resid   = b-res
        tempSol = resid./D + x
        diff    = sum(resid.^2)
        println(iter," ",diff/n)
    end
    return x
end

type TermStrVal
    str::AbstractString
    value::Float64
end

type TermLvlVal
    level::AbstractString
    value::Float64
end

type ModelTerm 
    iModel::Int64
    trmStr::AbstractString
    nFactors::Int64
    factors::Array{Symbol,1}
    str::Array{AbstractString,1}            # used to store the data for this term as strings
    val::Array{Float64,1}
    startPos::Int64                         # start pos in HMME
    nLevels::Int64                           
    X::SparseMatrixCSC{Float64,Int64}
    names::Array{Any,1}
end

type ResVar
    R0::Array{Float64,2}
    RiDict::Dict{BitArray{1},Array{Float64,2}}
end   

type MarkerMatrix
    X::Array{Float64,2}
    xArray::Array{Array{Float64,1},1}
    markerMeans::Array{Float64,2}
    mean2pq::Float64
    G::Array{Float64,1}
    centered::Bool
    function MarkerMatrix(X::Array{Float64,2},G::Array{Float64,1})
        markerMeans = center!(X) #centering
        p           = markerMeans/2.0
        mean2pq     = (2*p*(1-p)')[1,1]
        xArray      = get_column_ref(X)
        new(X,xArray,markerMeans,mean2pq,G,true)
    end
end

type MME
    modelVec::Array{AbstractString,1}
    modelTerms::Array{ModelTerm,1}
    modelTermDict::Dict{AbstractString,ModelTerm}
    lhsVec::Array{Symbol,1}
    covVec::Array{Symbol,1}
    pedTrmVec::Array{AbstractString,1}
    X
    ySparse
    mmeLhs
    mmeRhs
    ped
    Gi::Array{Float64,2}
    R::Array{Float64,2}
    Ai
    mmePos::Int64
    missingPattern
    resVar
    M
end

function mkDict(a)
    aUnique = unique(a)
    d = Dict()
    names = Array(Any,size(aUnique,1))
    for (i,s) in enumerate(aUnique)
        names[i] = s
        d[s] = i
    end
    return d,names
end

function getRi(resVar::ResVar,sel::BitArray{1})
    if haskey(resVar.RiDict,sel)
        return resVar.RiDict[sel]
    end
    n = size(resVar.R0,1)
    RZ = zeros(n,n)
    RZ[sel,sel] = inv(resVar.R0[sel,sel])
    resVar.RiDict[sel] = copy(RZ)
    return RZ
end

function getTerm(trmStr,m)
    trm = ModelTerm(m,string(m)*":"*trmStr,0,[],[],[],0,0,spzeros(0,0),[])
    factorVec = split(trmStr,"*")
    trm.nFactors = length(factorVec)
    trm.factors = [symbol(strip(f)) for f in factorVec]
    return trm
end

function initMME(models::AbstractString,R::Array{Float64,2})
    # returns an MME object for muilding the mme corresponding 
    # to the input string
    if models==""
        println("modelEquation is empty\n")
        return
    end
    modelVec = split(models,[';','\n'],keep=false)
    nModels  = size(modelVec,1)
    lhsVec   = Symbol[]
    modelTerms = ModelTerm[]
    dict = Dict{AbstractString,ModelTerm}()
    for (m,model) = enumerate(modelVec)
        lhsRhs = split(model,"=")
        lhsVec = [lhsVec;symbol(strip(lhsRhs[1]))]
        rhs = strip(lhsRhs[2])
        rhsVec = split(rhs,"+")    
        mTrms = [getTerm(strip(trmStr),m) for trmStr in rhsVec]
        modelTerms = [modelTerms; mTrms]
        for (i,trm) = enumerate(modelTerms) 
            dict[trm.trmStr] = modelTerms[i]
        end 
    end
    return MME(modelVec,modelTerms,dict,lhsVec,[],[],0,0,0,0,0,Array(Float64,1,1),R,0,1,0,0,0)
end 

function getData(trm::ModelTerm,df::DataFrame,mme::MME)
    nObs = size(df,1)
    trm.str = Array(AbstractString,nObs)
    trm.val = Array(Float64,nObs)
    if(trm.factors[1] == :intercept)                     # from Melanie's HW
        str = fill(string(trm.factors[1]),nObs)
        val = fill(1.0,nObs)
    else
        myDf = df[trm.factors]
        if trm.factors[1] in mme.covVec
            str = fill(string(trm.factors[1]),nObs)
            val = df[trm.factors[1]]
        else
            str = [string(i) for i in df[trm.factors[1]]]
            val = fill(1.0,nObs)
        end
        for i=2:trm.nFactors
            if trm.factors[i] in mme.covVec
                str = str .* fill("*"*string(trm.factors[i]),nObs)
                val = val .* df[trm.factors[i]]
            else
                str = str .* fill("*",nObs) .* [string(j) for j in df[trm.factors[i]]]
                val = val .* fill(1.0,nObs)
            end
        end
    end
    trm.str = str
    trm.val = val
end

getFactor1(str) = [strip(i) for i in split(str,"*")][1]

function getX(trm,mme::MME)
    pedSize = 0
    nObs  = size(trm.str,1)
    if trm.trmStr in mme.pedTrmVec
        trm.names   = PedModule.getIDs(mme.ped)
        trm.nLevels = length(mme.ped.idMap)
        xj = round(Int64,[mme.ped.idMap[getFactor1(i)].seqID for i in trm.str])
    else
        dict,trm.names  = mkDict(trm.str)
        trm.nLevels     = length(dict)
        xj    = round(Int64,[dict[i] for i in trm.str]) 
    end
    xi    = (trm.iModel-1)*nObs + collect(1:nObs)
    xv    = trm.val
    if mme.ped!=0
        pedSize = length(mme.ped.idMap)
        if trm.trmStr in mme.pedTrmVec
            # This is to ensure the X matrix for 
            # additive effect has the correct number of columns
            ii = 1         # adding a zero to
            jj = pedSize   # the last column in row 1
            vv = [0.0]
            xi = [xi;ii]
            xj = [xj;jj]
            xv = [xv;vv]
        end
    end 
    #make sure X has nObs*nModels rows
    nModels = size(mme.lhsVec,1)
    xi = [xi;1;nObs*nModels]
    xj = [xj;1;1]
    xv = [xv;0;0]
    trm.X = sparse(xi,xj,xv)
    trm.startPos = mme.mmePos
    mme.mmePos  += trm.nLevels
end

function mkRi(mme::MME,df::DataFrame)
    resVar = ResVar(mme.R,Dict())
    tstMsng = !isna(df[mme.lhsVec[1]])
    for i=2:size(mme.lhsVec,1)
        tstMsng = [tstMsng !isna(df[mme.lhsVec[i]])]
    end
    mme.missingPattern = tstMsng
    n    = size(tstMsng,2)
    nObs = size(tstMsng,1)
    ii = Array(Int64,nObs*n^2)
    jj = Array(Int64,nObs*n^2)
    vv = Array(Float64,nObs*n^2)
    pos = 1
    for i=1:size(tstMsng,1)
        sel = reshape(tstMsng[i,:],n)
        Ri  = getRi(resVar,sel)
        for ti=1:n
            tii = (ti-1)*nObs + i
            for tj=1:n
                tjj = (tj-1)*nObs + i
                ii[pos] = tii
                jj[pos] = tjj
                vv[pos] = Ri[ti,tj]
                pos += 1
            end
        end         
    end
    mme.resVar = resVar
    return sparse(ii,jj,vv)
end

function getMME(mme::MME, df::DataFrame)
    mme.mmePos = 1
    for trm in mme.modelTerms
        getData(trm,df,mme)
        getX(trm,mme)
    end
    n   = size(mme.modelTerms,1)
    trm = mme.modelTerms[1]
    X   = trm.X
    for i=2:n
        trm = mme.modelTerms[i]
        X = [X trm.X]
    end
    y = convert(Array,df[mme.lhsVec[1]],0.0)
    for i=2:size(mme.lhsVec,1)
        y    = [y; convert(Array,df[mme.lhsVec[i]],0.0)] 
    end
    N  = size(y,1)
    ii = 1:N
    jj = fill(1,N)
    vv = y
    ySparse = sparse(ii,jj,vv)
    nObs = size(df,1)
    Ri = mkRi(mme,df)
    mme.X = X
    mme.ySparse = ySparse 
    mme.mmeLhs = X'Ri*X
    mme.mmeRhs = X'Ri*ySparse
    if mme.ped != 0
        ii,jj,vv = PedModule.HAi(mme.ped)
        HAi = sparse(ii,jj,vv)
        mme.Ai = HAi'HAi
        addA(mme)
    end   
end

function getNames(mme)
    names = Array(String,0)
    for trm in mme.modelTerms
        for name in trm.names
            push!(names,trm.trmStr*": "*name)
        end
    end
    return names
end  

function covList(mme::MME, covStr::AbstractString)
    covVec = split(covStr," ",false) 
    mme.covVec = [symbol(i) for i in covVec]
    nothing
end

function getSolJ(mme::MME, df::DataFrame)
    if size(mme.mmeRhs)==() 
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    return [getNames(mme) Jacobi(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,0.3,tol=0.000001)]
end

function getSolG(mme::MME, df::DataFrame)
    if size(mme.mmeRhs)==() 
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    return [getNames(mme) GaussSeidel(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,tol=0.000001)]
end

function setAsRandom(mme::MME,randomStr::AbstractString,ped::PedModule.Pedigree, G::Array{Float64,2})
    pedTrmVec = split(randomStr," ",keep=false)
    res = []
    for trm in pedTrmVec
        for (m,model) = enumerate(mme.modelVec)
            strVec  = split(model,['=','+'])
            strpVec = [strip(i) for i in strVec]
            if trm in strpVec
                res = [res;string(m)*":"*trm]
            end
        end
    end
    mme.pedTrmVec = res
    mme.ped = ped
    mme.Gi = inv(G)
    nothing
end

function addA(mme::MME)
    pedTrmVec = mme.pedTrmVec
    for (i,trmi) = enumerate(pedTrmVec)
        pedTrmi  = mme.modelTermDict[trmi]
        startPosi  = pedTrmi.startPos
        endPosi    = startPosi + pedTrmi.nLevels - 1
        for (j,trmj) = enumerate(pedTrmVec)
            pedTrmj  = mme.modelTermDict[trmj]
            startPosj  = pedTrmj.startPos
            endPosj    = startPosj + pedTrmj.nLevels - 1 
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] = 
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] + mme.Ai*mme.Gi[i,j] 
        end
    end
end

function sampleMissingResiduals(mme,resVec)
    msngPtrn = mme.missingPattern
    n,k = size(msngPtrn)
    yIndex = collect(0:k-1)*n
    allTrue = fill(true,k)
    for i=1:n
        notMsng = reshape(msngPtrn[i,:,],k)
        if (notMsng!=allTrue)
            msng    = !notMsng
            nMsng   = sum(msng)
            resi    = resVec[yIndex+i][notMsng]
            Ri      = mme.resVar.RiDict[notMsng][notMsng,notMsng]
            Rc      = mme.R[msng,notMsng]
            L       = chol(mme.R[msng,msng] - Rc*Ri*Rc')'
            resVec[(yIndex+i)[msng]] = Rc*Ri*resi + L*randn(nMsng) 
            #resVec[yIndex+i][msng] = Rc*Ri*resi + L*randn(nMsng) this does not work!
        end
    end
end

ped = PedModule.mkPed("sim.ped");

dfGen = readtable("sim.gen", separator = ' ');

srand(1234)
Q = convert(Array{Float64,2},dfGen[:,collect(2:end)]);
α1 = randn(200)
α2 = randn(200)
a1 = Q*α1
a2 = Q*α2;

D = diagm(vec(sqrt(var([a1 a2],1))'))

R = cor([a1 a2])

G0 = D*R*D
G  = diag(G0)/200

R0 = diagm(vec(var([a1 a2],1)))
L  = chol(R0)
e  = L*randn(2,size(Q,1))
y = [a1 a2] + e'

df2 = DataFrame(Animal = dfGen[:,1], y1=y[:,1],y2=y[:,2])

head(df2)

models = "y1 = intercept ;
          y2 = intercept"
R = R0
mme = initMME(models,R)
nothing

function get_column(X,j)
    nrow,ncol = size(X)
    if j>ncol||j<0
        error("column number is wrong!")
    end
    indx = 1 + (j-1)*nrow
    ptr = pointer(X,indx)
    pointer_to_array(ptr,nrow)
end

function get_column_ref(X)
    ncol = size(X)[2]
    xArray = Array(Array{Float64,1},ncol)
    for i=1:ncol
        xArray[i] = get_column(X,i)
    end
    return xArray
end

function center!(X)
    nrow,ncol = size(X)
    colMeans = mean(X,1)
    BLAS.axpy!(-1,ones(nrow)*colMeans,X)
    return colMeans
end

function sample_effects_ycorr!(X,xArray,xpx,yCorr,α,meanAlpha,vRes,vEff,iIter) 
    nObs,nEffects = size(X)
    for j=1:nEffects
        x = xArray[j]
        rhs = dot(x,yCorr) + xpx[j]*α[j,1]
        lhs      = xpx[j] + λ
        invLhs   = 1.0/lhs
        mean     = invLhs*rhs
        oldAlpha = α[j,1]
        α[j]     = mean + randn()*sqrt(invLhs*vRes)
        BLAS.axpy!(oldAlpha-α[j,1],x,yCorr)
        meanAlpha[j] += (α[j] - meanAlpha[j])/iIter
    end
end


mme.M = MarkerMatrix(Q,G);

function sampleMCMC(nIter,mme,df;outFreq=100)
    getMME(mme,df)
    p = size(mme.mmeLhs,1)
    sol = fill(0.0,p)
    solMean = fill(0.0,p)
    GaussSeidel(mme.mmeLhs,sol,mme.mmeRhs,tol=0.000001) 
    ν = 10
    nObs    = size(df,1)
    nTraits = size(mme.lhsVec,1)
    νR0 = ν + nTraits
    R0 = mme.R
    PRes = R0*(νR0 - nTraits - 1)
    SRes   = zeros(Float64,nTraits,nTraits)
    R0Mean = zeros(Float64,nTraits,nTraits)
    if mme.ped != 0
        pedTrmVec = mme.pedTrmVec
        k = size(pedTrmVec,1)        
        νG0 = ν + k
        G0 = inv(mme.Gi)
        P = G0*(νG0 - k - 1)
        S = zeros(Float64,k,k)
        G0Mean = zeros(Float64,k,k)
    end

    for iter=1:nIter
        if iter%outFreq==0
            println("at sample: ",iter,"\n")
        end
        Gibbs(mme.mmeLhs,sol,mme.mmeRhs)
        # can make this more efficient by taking advantage of symmetry
        #if mme.ped != 0
            for (i,trmi) = enumerate(pedTrmVec)    
                pedTrmi  = mme.modelTermDict[trmi]
                startPosi  = pedTrmi.startPos
                endPosi    = startPosi + pedTrmi.nLevels - 1
                for (j,trmj) = enumerate(pedTrmVec)
                    pedTrmj  = mme.modelTermDict[trmj]
                    startPosj  = pedTrmj.startPos
                    endPosj    = startPosj + pedTrmj.nLevels - 1
                    S[i,j] = (sol[startPosi:endPosi]'*mme.Ai*sol[startPosj:endPosj])[1,1]
                end
            end
        #end
        resVec = mme.ySparse - mme.X*sol
        sampleMissingResiduals(mme,resVec)
        for traiti = 1:nTraits
            startPosi = (traiti-1)*nObs + 1
            endPosi   = startPosi + nObs - 1
            for traitj = traiti:nTraits
                startPosj = (traitj-1)*nObs + 1
                endPosj   = startPosj + nObs - 1
                SRes[traiti,traitj] = (resVec[startPosi:endPosi]'resVec[startPosj:endPosj])[1,1] 
                SRes[traiti,traitj] = SRes[traitj,traiti]
            end
        end
        R0 = rand(InverseWishart(νR0 + nObs, PRes + SRes))
        mme.R = R0
        Ri = mkRi(mme,df)
        X = mme.X
        mme.mmeLhs = X'Ri*X
        mme.mmeRhs = X'Ri*mme.ySparse
        if mme.ped != 0
            pedTrm1 = mme.modelTermDict[pedTrmVec[1]]
            q = pedTrm1.nLevels
            G0 = rand(InverseWishart(νG0 + q, P + S))
            mme.Gi = inv(G0)
            addA(mme)
            G0Mean  += (G0  - G0Mean )/iter
        end
        solMean += (sol - solMean)/iter
        R0Mean  += (R0  - R0Mean )/iter
    end
    output = Dict()
    output["posteriorMeanLocationParms"] = solMean
    #output["posteriorMeanG0"] = G0Mean
    output["posteriorMeanR0"] = R0Mean
    return output
end

res = sampleMCMC(200,mme,df2)
#res["posteriorMeanG0"]

mme.ped
