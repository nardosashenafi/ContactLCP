# using ContactLCP
using LaTeXStrings
using PyPlot

include("dynamics.jl")
include("../../src/lcp.jl")
include("../../src/solver.jl")

param = Float32[30.0, 5.0]
# param = Float32[0.0, 0.0]

sys  = RimlessWheel(Float32)
x0 = initialState(sys, 0.0f0, -0.5f0, 0.0f0, 0.0f0)
# x0 = Float32(initialState(sys, rand(-sys.α+0.1:0.05:0.0), rand(-2.0:0.05:-0.5), 0.0f0, 0.0f0))
# x0 = Float32([0.0, 0.26, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0])
lcp  = Lcp(Float32, sys)

X, t, Λn, Λt = fulltimestep(lcp, x0, param; Δt = 0.001f0, totalTimeStep = 5000)
plots(sys, X, t, Λn, Λt)

