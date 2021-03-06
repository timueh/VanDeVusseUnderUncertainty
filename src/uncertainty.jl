export 	quad_beta,
		buildBetaMixtureDensity,
		constructBasis,
		sampleFromBiMixture,
		verifyBasis

"""
	quad_beta(N::Int, α::Real, β::Real, w::Real)
returns ``N``-point quadrature for Beta distribution with shape parameters α, β, multiplied times w. 
"""
function quad_beta(N::Int, α::Real, β::Real, w::Real)
    op = Beta01OrthoPoly(N, α, β; addQuadrature=true)
    op.quad.nodes, w*op.quad.weights
end

"""
	buildBetaMixtureDensity(α::Vector, β::Vector, w::Vector)
returns density of beta mixture with vector of shape parameters α, β and weights w.
Note that the weights must sum to one --> the mixture density is a convex combination of standard Beta distribution densities.
"""
function buildBetaMixtureDensity(α::Vector, β::Vector, w::Vector)
	!(length(α) == length(β) == length(w)) && throw(DomainError((length(α), length(β), length(w)), "inconsistent lengths"))
	sum(w) != 1 && throw(DomainError(sum(w),"weights do not sum to one"))
	minimum(w) < 0 && throw(DomainError(minimum(w),"weights must be non-negative"))
    t -> sum( w_*PolyChaos.w_beta(t,α_,β_) for (α_, β_, w_) in zip(α, β, w) )
end

function verifyBasis(opq::AbstractOrthoPoly)
	deg, meas, ρ = opq.deg, opq.measure, opq.measure.w
	T2 = Tensor(2, opq)
	T2num_1 = [ T2.get([i,j]) for i in 0:deg, j in 0:deg]
	@show T2num_2 = [quadgk(x -> evaluate(i,x,opq)*evaluate(j,x,opq)*ρ(x),meas.dom[1],meas.dom[2])[1] for i in 0:deg, j in 0:deg ]
	T2num_2 - T2num_1
end

"""
	constructBasis(deg::Int, α::Vector, β::Vector, w::Vector)
construct orthogonal basis for Beta mixture model
"""
function constructBasis(deg::Int, α::Vector, β::Vector, w::Vector; Nrec::Int=2*deg)
	measure = Measure("MyMixture_Measure", buildBetaMixtureDensity(α,β,w), (0,1), false, Dict(:α=>α,:β=>β,:w=>w))
	a, b = mcdiscretization(Nrec, [n -> quad_beta(n, α[1], β[1], w[1]); n -> quad_beta(n, α[2], β[2], w[2])], gaussquad=true, Nmax=10000, ε=1e-4)
	OrthoPoly("MyMixtureOrthoPoly", deg, a, b, measure)
end

"""
	function sampleFromBiMixture(N::Int, α::Vector, β::Vector, w::Vector)
drawn ``N`` samples from a mixture of two Beta distributions
"""
function sampleFromBiMixture(N::Int, α::Vector, β::Vector, w::Vector)
	# check inputs
    length(w) != 2 && throw(DomainError(length(w), "method works for bi-variate mixtures only"))
    !(length(α) == length(β) == length(w)) && throw(DomainError((length(α), length(β), length(w)), "inconsistent lengths"))
    sum(w) != 1 && throw(DomainError(sum(w),"weights do not sum to one"))
    minimum(w) < 0 && throw(DomainError(minimum(w),"weights must be non-negative"))
    # begin sampling
    betas = [ Beta(a,b)  for (a,b) in zip(α,β) ]
    betaSamples = [ rand(beta, N) for beta in betas ]
    uniSamples = rand(N)
    mixtureSamples = Float64[]
    for n in 1:N
        ind = rand() < first(w) ? 1 : 2
        push!(mixtureSamples, betaSamples[ind][n])
    end
    mixtureSamples
end