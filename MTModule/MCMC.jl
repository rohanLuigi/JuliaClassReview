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
            println("at sample: ",iter)
            println(G0Mean)
        end
        Gibbs(mme.mmeLhs,sol,mme.mmeRhs)
        # can make this more efficient by taking advantage of symmetry
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
        resVec = mme.ySparse - mme.X*sol
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
        end
        solMean += (sol - solMean)/iter
        G0Mean  += (G0  - G0Mean )/iter
        R0Mean  += (R0  - R0Mean )/iter
    end
    output = Dict()
    output["posteriorMeanLocationParms"] = solMean
    output["posteriorMeanG0"] = G0Mean
    output["posteriorMeanR0"] = R0Mean
    return output
end
