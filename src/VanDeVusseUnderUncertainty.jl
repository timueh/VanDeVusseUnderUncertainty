__precompile__()

# unfortunately, compilation times are rather slow because of DifferentialEquations

module VanDeVusseUnderUncertainty
using 	JuMP, PolyChaos, LinearAlgebra, Ipopt,
		Distributions, DifferentialEquations, QuadGK

include("types.jl")
include("auxfuns.jl")
include("uncertainty.jl")
include("optimization.jl")

end # module
