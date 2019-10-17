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

Nsamples = 20
samples, Φ = generateSamples(Nsamples, unc)

A11_samples = vec(A[1,1,:]'*Φ)
x0_samples = x0*Φ

A_realizations = [ [A11 A[1,2,1]; A[2,1,1] A[2,2,1]] for A11 in A11_samples ]
B_realizations = [ B[:,:,1] for i in 1:Nsamples ]
x0_realizations = [ Vector(col) for col in eachcol(x0_samples) ]

sols = solveSystemForRealizations(A_realizations, B_realizations, x0_realizations, usol, N)

##
using PyPlot, LaTeXStrings
t_x, t_u = 0:δ:T, 0:δ:T-δ
close("all")
subplot(2,2,1)
grid(true)
fill_between(t_x, μ1 + 3*σ1, μ1 - 3*σ1, alpha=0.1, color="k", edgecolor="k")
plot(t_x, μ1)
[ plot(t_x, x1sol, "--") for (x1sol, x2sol) in sols]
xlabel(L"t"); ylabel(L"x_1(t)")
ylim([0, 0.55])
subplot(2,2,2)
grid(true)
xlabel(L"t"); ylabel(L"x_2(t)")
fill_between(t_x, μ2 + 3*σ2, μ2 - 3*σ2, alpha=0.1, color="k", edgecolor="k")
plot(t_x, μ2)
[ plot(t_x, x2sol, "--") for (x1sol, x2sol) in sols]
plot(t_x, x2max*ones(size(t_x)), "--r")
ylim([0, 0.20])
subplot(2,2,3)
grid(true)
xlabel(L"t"); ylabel(L"u(t)")
plot(t_u, usol[1,:])




