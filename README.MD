# Optimal control under uncertainty

*This package is work in progress*

Given a discrete-time linear time-invariant system subject to parametric uncertainties, this code uses polynomial chaos expansion to compute optimal input trajectories such that the expected quadratic cost is minimized subject to individual chance constraints.

The example is a linearized van-de-Vusse reactor.

We consider three distinct sources of uncertainty:
    - the reaction rate is modelled as a Beta mixture model, and
    - the initial condition for each concentratio is uncertain.