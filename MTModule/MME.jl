module MTModule

using PedModule
using DataFrames

include("solver.jl")
include("types.jl")
include("functions.jl")
include("MCMC.jl")

export initMME
export setAsRandom
export sampleMCMC

end
