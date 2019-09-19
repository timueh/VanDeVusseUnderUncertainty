__precompile__()

module OptimalControlUnderUncertainty
using JuMP, PolyChaos, LinearAlgebra, Ipopt, Distributions

include("types.jl")
include("auxfuns.jl")
include("uncertainty.jl")
include("optimization.jl")

end # module
