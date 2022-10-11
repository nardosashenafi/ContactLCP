# using ContactLCP
using LaTeXStrings
using PyPlot

include("dynamics.jl")
include("../../src/lcp.jl")
include("../../src/solver.jl")

param = Float64[20.0, 5.0]
# param = Float64[0.0, 0.0]

sys  = RimlessWheel(Float64)
# x0 = initialState(sys, 0.0, -0.5f0, 0.0f0, 0.0f0)
x0 =   [0.08035215340930911
0.2482699308407366
0.29506445059426006
-0.31420206637617604
-0.010289938728842103
0.0033441274048094227
0.029285686470857678
0.041614251929509606]
# x0 = [0.0, 0.26, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0]
lcp  = Lcp(Float64, sys)


X, t, Λn, Λt = fulltimestep(lcp, x0, param; Δt = 0.001, totalTimeStep = 5000)
plots(sys, X, t, Λn, Λt)

