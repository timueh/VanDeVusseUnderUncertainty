using JuMP, PolyChaos, LinearAlgebra, Ipopt, Distributions



function getSolution(x::Matrix, op::AbstractOrthoPoly)
	mu  = [ mean(Vector(row), op) for row in eachrow(x) ]
	sig = [ std(Vector(row), op) for row in eachrow(x) ]
	mu, sig
end

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



nx, nu = 2, 1
d = 4
α, β, w = [2, 4], [4.5, 1.5], [0.3, 0.7]

lb, ub = 0.923, 0.963

op_k = constructBasis(d,α,β,w)
op_x1 = GaussOrthoPoly(d; Nrec=2*d)
op_x2 = GaussOrthoPoly(d; Nrec=2*d)
mop = MultiOrthoPoly([op_k, op_x1, op_x2], d)

L = dim(mop)
T2, T3 = Tensor(2, mop), Tensor(3, mop)

function tensor(i,j,k)
	T3.get([i-1, j-1, k-1]) / T2.get([k-1, k-1])
end

A = zeros(nx, nx, L)
A[:,:,1] = [ 0 0; 0.088 0.819 ]
A[1,1,1:2] = convert2affinePCE(lb, ub-lb, first(op_k.α))

B = zeros(nx, nu, L)
B[:,:,1] = [ -0.005; -0.002 ]


δ = 0.002
T = Int(0.2/δ)

x0 = zeros(nx, L)
x0[1,[1,3]] = convert2affinePCE(0.5, 0.1/3*0.5, op_x1) 
x0[2,[1,4]] = convert2affinePCE(0.1, 0.1/3*0.1, op_x2)
Q, R = diagm(0 => [1; 1]), diagm(0 => 1e0*ones(nu))
x2max, λ = 0.17, 1.618
umax, Δumax = 1.0, 0.2

model = Model(with_optimizer(Ipopt.Optimizer))
@variable(model, x[1:nx, 1:T+1, 1:L])
@variable(model, var2[1:T+1] >= 0)
@variable(model, u[1:nu, 1:T])
@constraint(model, Dynamics[t in 1:T,k in 1:L], x[:,t+1,k] .== sum( A[:,:,l]*x[:,t,m]*tensor(l,m,k) for l in 1:L, m in 1:L) + B[:,:,k]*u[:,t])
@constraint(model, InitialCondition[k in 1:L], x[:,1,k] .== x0[:,k])
@constraint(model, Variance[t in 1:T+1], var2[t] == sum(x[2,t,l]^2*T2.get([l-1, l-1]) for l in 2:L))
@NLconstraint(model, ChanceConstraint[t in 1:T+1], λ^2*var2[t] <= (x2max - x[2,t,1])^2 )
# @constraint(model, Actuator[t in 2:T], -Δumax .<= u[:,t] - u[:,t-1] .<= Δumax)
@objective(model, Min, sum(dot(x[:,t,m], Q*x[:,t,l])*T3.get([m-1,l-1,0]) for t in 1:T+1, m in 1:L, l in 1:L) + sum(dot(u[:,t], R*u[:,t]) for t in 1:T) )
optimize!(model)

x1sol, x2sol = value.(x)[1,:,:], value.(x)[2,:,:]

μ1, σ1 = getSolution(x1sol, mop)
μ2, σ2 = getSolution(x2sol, mop)

##
using PyPlot, LaTeXStrings
t_x, t_u = δ*(0:T), δ*(0:T-1)

close("all")
subplot(2,2,1)
grid(true)
fill_between(t_x, μ1 + 3*σ1, μ1 - 3*σ1, alpha=0.1, color="k", edgecolor="k")
plot(t_x, μ1)
xlabel(L"t"); ylabel(L"x_1(t)")
subplot(2,2,2)
grid(true)
xlabel(L"t"); ylabel(L"x_2(t)")
fill_between(t_x, μ2 + 3*σ2, μ2 - 3*σ2, alpha=0.1, color="k", edgecolor="k")
plot(t_x, μ2)
plot(t_x, x2max*ones(size(t_x)), "--r")
subplot(2,2,3)
grid(true)
xlabel(L"t"); ylabel(L"u(t)")
plot(t_u, value.(u)[1,:])


