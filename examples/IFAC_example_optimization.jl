using 	VanDeVusseUnderUncertainty, JuMP, PolyChaos, LinearAlgebra,
		Ipopt, Distributions, QuadGK


degree = 4
α, β, w = [2, 4], [4.5, 1.5], [0.3, 0.7]

op_k = constructBasis(degree,α,β,w)
op_x1 = GaussOrthoPoly(degree; Nrec=2*degree)
op_x2 = GaussOrthoPoly(degree; Nrec=2*degree)
mop = MultiOrthoPoly([op_k, op_x1, op_x2], degree)
unc = Uncertainty(mop)
L = VanDeVusseUnderUncertainty.dim(unc)
lb, ub = 0.923, 0.963

δ, T = 0.002, 0.2
N = Int(T/δ)

nx, nu = 2, 1
A = zeros(nx, nx, L)
A[:,:,1] = [ 0 0; 0.088 0.819 ]
A[1,1,1:2] = convert2affinePCE(lb, ub-lb, first(op_k.α))

B = zeros(nx, nu, L)
B[:,:,1] = [ -0.005; -0.002 ]

x0 = zeros(nx, L)
x0[1,[1,3]] = convert2affinePCE(0.5, 0.1/3*0.5, op_x1)
x0[2,[1,4]] = convert2affinePCE(0.1, 0.1/3*0.1, op_x2)

sys = LTISystem(A,B,x0,δ)

cost = Cost(diagm(0 => [1; 1]), diagm(0 => 1e0*ones(nu)))

x2max, λ = 0.17, 1.618
model = createOptimizationProblem(sys, cost, unc, N, x2max=x2max, λ=λ)
optimize!(model)
μ1, σ1, μ2, σ2, usol = getSolution(model, unc)


##
include("optiPlot.jl")
