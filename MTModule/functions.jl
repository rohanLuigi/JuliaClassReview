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
    return MME(modelVec,modelTerms,dict,lhsVec,[],[],0,0,0,0,0,Array(Float64,1,1),R,0,1)
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

function mkRi(mme::MME,df::DataFrame)
    resVar = ResVar(mme.R,Dict())
    tstMsng = !isna(df[mme.lhsVec[1]])
    for i=2:size(mme.lhsVec,1)
        tstMsng = [tstMsng !isna(df[mme.lhsVec[i]])]
    end
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
    mme.mmeLhs = X'Ri*X
    mme.mmeRhs = X'Ri*ySparse
    if mme.ped != 0
        ii,jj,vv = PedModule.HAi(mme.ped)
        HAi = sparse(ii,jj,vv)
        mme.Ai = HAi'HAi
        addA(mme::MME)
    end
end

function getNames(mme)
    names = Array(AbstractString,0)
    for trm in mme.modelTerms
        for name in trm.names
            push!(names,trm.trmStr*": "*name)
        end
    end
    return names
end

function covList(mme::MME, covStr::AbstractString)
    covVec = split(covStr," ",keep=false)
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
