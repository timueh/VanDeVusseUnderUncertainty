using VanDeVusseUnderUncertainty, DifferentialEquations, PolyChaos, PyPlot, LaTeXStrings

function ODEdeterministic(du, u, p, t)
	k_rate, input = p[1], p[2]
    du[1] = -k_rate[1]*u[1] - k_rate[3]*u[1]^2 - u[1]*input
    du[2] =  k_rate[1]*u[1] - k_rate[2]*u[2]   - u[2]*input
end

function ODEgalerkin(du, u, p, t)
    x1, x2 = u[1:L], u[L+1:end]
    du[1:L] = [ -x1[k]*p - sum((k1[l]*x1[m] + k_rate[3]*x1[l]*x1[m])*T3.get([l-1, m-1, k-1])/T2.get([k-1, k-1]) for l in inds for m in inds) for k in inds]
    du[L+1:end] = [ -x2[k]*p + sum((k1[l]*x1[m]-k2[l]*x2[m])*T3.get([l-1, m-1, k-1])/T2.get([k-1, k-1]) for l in inds, m in inds) for k in inds ]
end

# find out reference for those numbers --> look up references in Paulson, Streif, Mesbah. "Stability for Receding-horizon Stochastic Model Predictive Control"
# reaction parameters
k0 = 1e9*[12.87; 12.87; 9.043]
T0 = 104.9
E = [9758.3; 9758.3; 8560]
k_rate = @. k0 * exp(-E/(T0 + 273.15))
k_rate = [50; 100; 10 ] # Constrained Linear Quadratic Regulation, Pierre O. M. Scokaert and James B. Rawlings
# initial condition
c0 = [0.5; 0.1]
tend, Δt = 0.1, 0.001


prob = ODEProblem(ODEdeterministic, c0, (0, tend), (k_rate, 0.1))
sol = DifferentialEquations.solve(prob; saveat=0:Δt:tend);

# uncertainty propagation
degree = 4
opq_1, opq_2 = Uniform01OrthoPoly(degree, Nrec = 4*degree), Uniform01OrthoPoly(degree, Nrec = 4*degree)
mop = MultiOrthoPoly([opq_1, opq_2], degree)
T2, T3 = Tensor(2, mop), Tensor(3, mop)
L = PolyChaos.dim(mop)
inds = 1:L
k1, k2 = zeros(L), zeros(L)
k1[1:2], k2[[1,3]] = convert2affinePCE(k_rate[1], 0.1*k_rate[1], opq_1, kind="μσ") , convert2affinePCE(k_rate[2], 0.1*k_rate[2], opq_2, kind="μσ")

e = zeros(L)
e[1] = 1

probGalerkin = ODEProblem(ODEgalerkin, kron(c0,e), (0, tend), 0.1)
sol = DifferentialEquations.solve(probGalerkin; saveat=0:Δt:tend)

t_x = sol.t

μ1, σ1, μ2, σ2 = getSolution(sol, mop)

# compare against sampled solution
N_samples = 20
k1_samples = samplePCE(N_samples, k1, mop)
k2_samples = samplePCE(N_samples, k2, mop)


function ODEProblemFunction(k1,k2)
	ODEProblem(ODEdeterministic,c0,(0,tend),([k1, k2, k_rate[3]], 0.1))
end

deterministic_solutions = solveSystemForRealizations(ODEProblemFunction, k1_samples, k2_samples, 0:Δt:tend)


width, height = 3.6, 1.8 # in inches!
close("all")
PyPlot.rc("text", usetex=true)
PyPlot.rc("font", size=9)
fig1 = figure(1, frameon=true, tight_layout=true, figsize=(width, height))
# fig1.tight_layout()
subplot(1,2,1)
PyPlot.grid(true)
fill_between(t_x, μ1 + 3*σ1, μ1 - 3*σ1, alpha=0.1, color="k", edgecolor="k")
[ plot(t_x, x1sol, "--",linewidth=.7) for (x1sol, x2sol) in deterministic_solutions ]
plot(t_x, μ1)
xlabel(L"t"); ylabel(L"c_{A}(t)")
#savefig("/home/ws/oz8348/ifac_polychaospaper/propagation_x1.pdf", format="pdf")

subplot(1,2,2)
PyPlot.grid(true)
fig2.tight_layout()
fill_between(t_x, μ2 + 3*σ2, μ2 - 3*σ2, alpha=0.1, color="k", edgecolor="k")
[ plot(t_x, x2sol, "--",linewidth=.7) for (x1sol, x2sol) in deterministic_solutions ]
plot(t_x, μ2)
xlabel(L"t"); ylabel(L"c_{B}(t)")
savefig("/home/ws/oz8348/ifac_polychaospaper/propagation.pdf", format="pdf")
