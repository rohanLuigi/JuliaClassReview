using Distributions

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

function getNames(mme)
    names = Array(AbstractString,0)
    for trm in mme.modelTerms
        for name in trm.names
            push!(names,trm.trmStr*": "*name)
        end
    end
    return names
end

function getTerm(trmStr) #only declaration of trm is changed
    trm = ModelTerm(trmStr,0,[],[],[],0,0,spzeros(0,0),[])
    factorVec = split(trmStr,"*")
    trm.nFactors = length(factorVec)
    trm.factors = [symbol(strip(f)) for f in factorVec]
    return trm
end

function initMME(modelEquation::AbstractString,R::Float64) #Difference: made modelTerms into a dictionary
    if modelEquation==""
        error("modelEquation is empty\n")
    end
    lhsRhs = split(modelEquation,"=")
    lhs = symbol(strip(lhsRhs[1]))
    rhs = strip(lhsRhs[2])
    rhsVec = split(rhs,"+")
    dict = Dict{AbstractString,ModelTerm}()
    modelTerms = [getTerm(strip(trmStr)) for trmStr in rhsVec]
    for (i,trm) = enumerate(modelTerms)
        dict[trm.trmStr] = modelTerms[i]
    end
    return MME(modelEquation,modelTerms,dict,lhs,R)
end

function getData(trm,df::DataFrame,mme::MME)
    nObs = size(df,1)
    trm.str = Array(AbstractString,nObs)
    trm.val = Array(Float64,nObs)

    if trm.factors[1] == :intercept ##modified from Melanie's HW ##new
        str = fill(string(trm.factors[1]),nObs)
        val = fill(1.0,nObs)
        trm.str = str
        trm.val = val
        return
    end

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
            str = str .* fill(" x "*string(trm.factors[i]),nObs)
            val = val .* df[trm.factors[i]]
        else
            str = str .* fill(" x ",nObs) .* [string(j) for j in df[trm.factors[i]]]
            val = val .* fill(1.0,nObs)
        end
    end
    trm.str = str
    trm.val = val
end

function covList(mme::MME, covStr::AbstractString)
    covVec = split(covStr," ",keep=false)
    mme.covVec = [symbol(i) for i in covVec]
    nothing
end

getFactor1(str) = [strip(i) for i in split(str,"x")][1] #using in may be better. maybe age*animal

function getX(trm::ModelTerm,mme::MME)
    pedSize = 0
    nObs  = size(trm.str,1)
    if trm.trmStr in mme.pedTrmVec #"Animal"
        trm.names   = PedModule.getIDs(mme.ped)
        trm.nLevels = length(mme.ped.idMap)
        xj = round(Int64,[mme.ped.idMap[getFactor1(i)].seqID for i in trm.str]) #column index
    else
        dict,trm.names  = mkDict(trm.str)
        trm.nLevels     = length(dict)
        xj    =  round(Int64,[dict[i] for i in trm.str])
    end
    xi    = 1:nObs  #row index
    xv    = trm.val #value
    if mme.ped!=0 #Because some animals in pedigree may not in the data
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
    trm.X = sparse(xi,xj,xv)
    trm.startPos = mme.mmePos
    mme.mmePos  += trm.nLevels
end

function getMME(mme::MME, df::DataFrame)
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
    y    = df[mme.lhs]
    nObs = size(y,1)
    ii = 1:nObs
    jj = fill(1,nObs)
    vv = y
    nRowsX = size(X,1)
    if nRowsX > nObs  ###?????? nRowsX=nObs
        ii = [ii,nRowsX]
        jj = [jj,1]
        vv = [vv,0.0]
    end
    ySparse = sparse(ii,jj,vv)
    mme.X = X
    mme.ySparse = ySparse
    mme.mmeLhs = X'X
    mme.mmeRhs = X'ySparse
    if mme.ped != 0
        ii,jj,vv = PedModule.HAi(mme.ped)#cholesky??
        HAi = sparse(ii,jj,vv)
        mme.Ai = HAi'HAi
        addA(mme::MME)
    end
    addLambdas(mme)
end

function setAsRandom(mme::MME,randomStr::AbstractString,ped::PedModule.Pedigree, G::Array{Float64,2})
    mme.pedTrmVec = split(randomStr," ",keep=false)
    mme.ped = ped
    mme.GiOld = zeros(G)
    mme.GiNew = inv(G)
    nothing
end

function setAsRandom(mme::MME,randomStr::AbstractString, vc::Float64, df::Float64)
    trm  = mme.modelTermDict[randomStr]
    scale = vc*(df-2)/df
    randomEffect = RandomEffect(trm,0.0,vc,df,scale)
    push!(mme.rndTrmVec,randomEffect)
end

function addA(mme::MME)
    pedTrmVec = mme.pedTrmVec
    #pedTrm = mme.modelTermDict[mme.pedTrmVec[1]]
    for (i,trmi) = enumerate(pedTrmVec)
        pedTrmi  = mme.modelTermDict[trmi]
        startPosi  = pedTrmi.startPos
        endPosi    = startPosi + pedTrmi.nLevels - 1
        for (j,trmj) = enumerate(pedTrmVec)
            pedTrmj  = mme.modelTermDict[trmj]
            startPosj  = pedTrmj.startPos
            endPosj    = startPosj + pedTrmj.nLevels - 1
            
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] =
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] + mme.Ai*(mme.GiNew[i,j]*mme.RNew - mme.GiOld[i,j]*mme.ROld)
        end
    end
end

function addLambdas(mme::MME)
    for effect in  mme.rndTrmVec
        trmi       = effect.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        lambdaDiff = mme.RNew/effect.vcNew - mme.ROld/effect.vcOld
        mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] = 
        mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] + eye(trmi.nLevels)*lambdaDiff
    end   
end

function getSolJ(mme::MME, df::DataFrame)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    return [getNames(mme) Jacobi(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,0.3,tol=0.000001)]
end

function getSolG(mme::MME, df::DataFrame;outFreq=10)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    return [getNames(mme) GaussSeidel(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,tol=0.000001,output=outFreq)]
end

function getSolGibbs(mme::MME, df::DataFrame;nIter=50000,outFreq=100)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    return [getNames(mme) Gibbs(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,mme.R,nIter,outFreq=outFreq)]
end


function sampleMCMC(nIter,mme,df;outFreq=100)
    if size(mme.mmeRhs)==() 
        MMEModule.getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    sol = zeros(p)
    solMean = zeros(p)
    
    initSampleArrays(mme,nIter)
    
    vRes=mme.RNew
    nuRes=4
    scaleRes   = vRes*(nuRes-2)/nuRes 

    if mme.M!=0
        nObs,nLoci = size(mme.M.X)
        α  = zeros(Float64,nLoci)
        meanAlpha = zeros(Float64,nLoci)
        mArray = mme.M.xArray
        mpm = [dot(mme.M.X[:,i],mme.M.X[:,i]) for i=1:size(mme.M.X,2)]   
        M = mme.M.X
        dfEffectVar=4
        vEff=mme.M.G/mme.M.mean2pq 
        scaleVar   = vEff*(dfEffectVar-2)/dfEffectVar        # scale factor for locus effects
    end
    
    ycorr = vec(full(mme.ySparse))
    
    ν = 10
    if mme.ped != 0
        pedTrmVec = mme.pedTrmVec
        k = size(pedTrmVec,1)
        νG0 = ν + k
        G0 = inv(mme.GiNew)
        P = G0*(νG0 - k - 1)
        S = zeros(Float64,k,k)
        G0Mean = zeros(Float64,k,k)
    end

    
    for iter=1:nIter
        #sample non-marker part
        ycorr = ycorr + mme.X*sol
        rhs = mme.X'ycorr #

        MMEModule.Gibbs(mme.mmeLhs,sol,rhs,vRes)
        ycorr = ycorr - mme.X*sol

        solMean += (sol - solMean)/iter
        
        #sample marker
        if mme.M!=0
            MMEModule.sample_effects_ycorr!(M,mArray,mpm,ycorr,α,meanAlpha,vRes,vEff,iter)
        end
 
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
        
        if mme.ped != 0
            pedTrm1 = mme.modelTermDict[pedTrmVec[1]]
            q = pedTrm1.nLevels
            G0 = rand(InverseWishart(νG0 + q, P + S)) #ν+q?
            mme.GiOld = copy(mme.GiNew)
            mme.GiNew = inv(G0)
            MMEModule.addA(mme)
        end
       
        sampleVCs(mme)
        
        addLambdas(mme::MME)
        
        mme.ROld = mme.RNew
        mme.RNew = sampleVariance(ycorr, length(ycorr), nuRes, scaleRes)
        outputSamples(mme,sol,iter)
        if mme.M!=0
            vEff = sampleVariance(α, nLoci, dfEffectVar, scaleVar)
        end

        if iter%outFreq==0
            println("at sample: ",iter)
        end
    end
    output = Dict()
    output["posteriorMeanLocationParms"] = [MMEModule.getNames(mme) solMean]
    if mme.M!=0
        output["posteriorMeanMarkerEffects"] = meanAlpha
    end
    for i in  mme.outputSamplesVec
        trmi   = i.term
        trmStr = trmi.trmStr
        output["MCMCSamples:"*trmStr] = i.sampleArray
    end
    return output
end

function sampleVCs(mme::MME)
    for effect in  mme.rndTrmVec
        trmi       = effect.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        x          = sol[startPosi:endPosi]
        effect.vcOld  = effect.vcNew
        effect.vcNew  = sampleVariance(x,trmi.nLevels, effect.df, effect.scale)
    end   
end

function outputSamplesFor(mme::MME,trmStr::AbstractString)
    trm  = mme.modelTermDict[trmStr]
    samples = MCMCSamples(trm,Array(Float64,1,1))
    push!(mme.outputSamplesVec,samples)
end

function initSampleArrays(mme::MME,niter)
    for i in  mme.outputSamplesVec
        trmi = i.term
        i.sampleArray = zeros(niter,trmi.nLevels)
    end
end

function outputSamples(mme::MME,sol,iter::Int64)
    for i in  mme.outputSamplesVec
        trmi = i.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        i.sampleArray[iter,:] = sol[startPosi:endPosi]
    end
end


    
