# using ContactLCP
using LaTeXStrings
using PyPlot

include("dynamics.jl")
include("../../src/lcp.jl")
include("../../src/solver.jl")

param = Float64[100, 20]

sys  = RimlessWheel(Float64)
lcp  = Lcp(Float64, sys)


X, t, Λn, Λt = fulltimestep(lcp, sys.x0, param; Δt = 0.001, totalTimeStep = 100)
plots(X, t, Λn, Λt)

