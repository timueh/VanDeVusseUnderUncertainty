using PyPlot, LaTeXStrings



t_x, t_u = 0:δ:T, 0:δ:T-δ
close("all")
PyPlot.rc("font", size=9)
fig1 = figure(1, frameon=true, tight_layout=true, figsize=(3.6, 1.8))

#subplot(1,2,1)
#grid(true)
#fill_between(t_x, μ1 + 3*σ1, μ1 - 3*σ1, alpha=0.1, color="k", edgecolor="k")
#plot(t_x, μ1)
#[ plot(t_x, x1sol, "--") for (x1sol, x2sol) in sols]
#xlabel(L"t"); ylabel(L"x_1(t)")
subplot(1,2,1)
grid(true)
xlabel(L"t"); ylabel(L"x_2(t)")
fill_between(t_x, μ2 + 3*σ2, μ2 - 3*σ2, alpha=0.1, color="k", edgecolor="k")
plot(t_x, μ2)
[ plot(t_x, x2sol, "--") for (x1sol, x2sol) in sols]
plot(t_x, x2max*ones(size(t_x)), "--r")
subplot(1,2,2)
grid(true)
xlabel(L"t"); ylabel(L"u(t)")
plot(t_u, usol[1,:])
#subplots_adjust(wspace=.2,hspace=.2)
savefig("/home/ws/oz8348/ifac_polychaospaper/optimization.pdf", format="pdf")

meas_plot1 = Beta01Measure(α[1],β[1])
meas_plot2 = Beta01Measure(α[2],β[2])

fig2 = figure(2, frameon=true, tight_layout=true, figsize=(3.3, 1.9))
subplot(1,2,1)
xgrid = 0:0.001:1
plot(xgrid, op_k.measure.w.(xgrid))
grid(true)
xlabel(L"\tau"); ylabel(L"\rho(\tau)")

function build_myfun(n::Int)
    t -> evaluate(n, t, op_k)
end

tspan = 0:0.01:1
subplot(1,2,2)
grid(true)
xlabel(L"\tau"); ylabel(L"\phi_k(\tau)")
[ plot(tspan, build_myfun(d).(tspan)) for d in 1:4 ]
PyPlot.ylim(-0.1, 0.3)

savefig("/home/ws/oz8348/ifac_polychaospaper/optimizationBetaPDF.pdf", format="pdf")
