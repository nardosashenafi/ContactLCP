
# using ContactLCP
using LaTeXStrings
using PyPlot

include("bouncing_ball.jl")
include("../src/lcp.jl")
include("../src/solver.jl")

sys  = BouncingBall(Float64)
lcp  = Lcp(sys, Float64)

X, t, Λn, Λt = fulltimestep(lcp, sys.x0, Float64; totalTimeStep=1500)
plots(X, t, Λn, Λt)