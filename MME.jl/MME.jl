include("solver.jl")

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

type ModelTerm 
    trmStr::AbstractString
    nFactors::Int64
    factors::Array{Symbol,1}
    data::Array{AbstractString,1}
    X::SparseMatrixCSC{Float64,Int64}
    names::Array{Any,1}
end

function getTerm(trmStr)
    trm = ModelTerm(trmStr,0,[],[],spzeros(0,0),[])
    if length(trmStr)==1
        trm.nFactors = 1
        trm.factors  = [symbol(strip(trmStr))]
    else
        factorVec = split(trmStr,"*")
        trm.nFactors = length(factorVec)
        trm.factors = [symbol(strip(f)) for f in factorVec]
    end
    return trm
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