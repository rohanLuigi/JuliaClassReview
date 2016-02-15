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

type RandomEffect
    term::ModelTerm
    vcOld::Float64
    vcNew::Float64
    df::Float64
    scale::Float64
end

type MME
    modelEquation::AbstractString
    modelTerms::Array{ModelTerm,1}
    modelTermDict::Dict{AbstractString,ModelTerm}
    lhs::Symbol
    covVec::Array{Symbol,1}
    pedTrmVec::Array{AbstractString,1}
    rndTrmVec::Array{RandomEffect,1}
    X
    ySparse
    mmeLhs
    mmeRhs
    ped
    GiOld::Array{Float64,2}
    GiNew::Array{Float64,2}
    ROld
    RNew
    Ai
    mmePos::Int64
    M #marker genotypes
    function MME(modelEquation,modelTerms,dict,lhs,R)
        new(modelEquation,modelTerms,dict,lhs,[],[],[],0,0,0,0,0,Array(Float64,1,1),Array(Float64,1,1),0.0,R,0,1,0)
    end
end
