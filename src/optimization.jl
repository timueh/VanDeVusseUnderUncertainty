export 	createOptimizationProblem

function createOptimizationProblem(sys::LTISystem, cost::Cost, unc::Uncertainty, T::Int; x2max::Real=0.17, λ::Real=1.618)
	A, B, nx, nu = sys.A, sys.B, sys.nx, sys.nu
	x0 = sys.x0
	Q, R = cost.Q, cost.R

	T2, T3 = unc.T2, unc.T3
	tensor = _buildTensorFun(unc)
	L = OptimalControlUnderUncertainty.dim(unc)

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
	model
end