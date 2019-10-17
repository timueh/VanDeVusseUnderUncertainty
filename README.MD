# Optimal control under uncertainty

*This package is work in progress*

This package studies a Van de Vusse reactor under uncertainty. There are three examples:
- `IFAC_example_basisConstruction.jl` which constructs the basis for a Beta mixture model,
- `IFAC_example_propagation.jl` which propagates uncertainties through the system, and
- `IFAC_example_optimization.jl` which solves an optimal control problem under uncertainty.

The equations are based on [this paper](https://ieeexplore.ieee.org/abstract/document/7170854/).

All uncertainties are modeled as continuous random variables of finite variance. Polynomial chaos expansion is employed to propagate uncertainties. We make use of the Julia package [`PolyChaos.jl`](https://github.com/timueh/PolyChaos.jl).