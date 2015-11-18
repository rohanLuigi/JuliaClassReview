type ModelTerm
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

type MME
    modelEquation::AbstractString
    modelTerms::Array{ModelTerm,1}
    modelTermDict::Dict{AbstractString,ModelTerm}
    lhs::Symbol
    covVec::Array{Symbol,1}
    pedTrmVec::Array{AbstractString,1}
    X
    ySparse
    mmeLhs
    mmeRhs
    ped
    Gi::Array{Float64,2}
    R
    Ai
    mmePos::Int64
    M #marker genotypes
end
