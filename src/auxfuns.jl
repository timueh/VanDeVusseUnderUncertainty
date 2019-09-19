export 	dim,
		getSolution,
		solveSystem,
		generateSamples,
		solveSystemForRealizations

function _issquare(A::Matrix)::Bool
	m, n = size(A)
	m == n
end

function _buildTensorFun(T2::Tensor, T3::Tensor)
	(i,j,k) -> T3.get([i-1, j-1, k-1]) / T2.get([k-1, k-1])
end

_buildTensorFun(unc::Uncertainty) = _buildTensorFun(unc.T2, unc.T3)

dim(unc::Uncertainty) = PolyChaos.dim(unc.op)

function getSolution(x::Matrix, op::AbstractOrthoPoly)
	mu  = [ mean(Vector(row), op) for row in eachrow(x) ]
	sig = [ std(Vector(row), op) for row in eachrow(x) ]
	mu, sig
end

getSolution(x::Matrix, unc::Uncertainty) = getSolution(x, unc.op)

function getSolution(model::Model,unc::Uncertainty)
	x1sol, x2sol = value.(model[:x])[1,:,:], value.(model[:x])[2,:,:]
	usol = value.(model[:u])
	μ1, σ1 = getSolution(x1sol, unc)
	μ2, σ2 = getSolution(x2sol, unc)
	μ1, σ1, μ2, σ2, usol
end

function solveSystem(A::Matrix, B::Matrix, x0::Vector, u, T)
	x = [x0]
	for t in 1:T
		push!(x, vec(A*x[t] + B*u[t]))
	end
	[ x_[1] for x_ in x ], [ x_[2] for x_ in x ]
end

function solveSystemForRealizations(A::Vector, B::Vector, x0::Vector, u, N)
	sols = []
	for (A_,B_,x0_) in zip(A,B,x0)
		push!(sols, solveSystem(A_,B_,x0_,u,N))
	end
	sols
end

function generateSamples(N::Int, op::AbstractOrthoPoly)
	N <= 0 && throw(DomainError(N,"number of samples must be positive"))
	α, β, w = op.uni[1].measure.pars[:α], op.uni[1].measure.pars[:β], op.uni[1].measure.pars[:w]
	op_x1 = op.uni[2]
	op_x2 = op.uni[3]
	samples = [sampleFromBiMixture(N,α,β,w) sampleMeasure(N, op_x1) sampleMeasure(N, op_x2)]
	Φ = evaluate(samples, op)
	samples, Φ
end

generateSamples(N::Int, unc::Uncertainty) = generateSamples(N, unc.op)