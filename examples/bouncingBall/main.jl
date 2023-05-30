
using ContactLCP
using LaTeXStrings
using PyPlot

include("bouncing_ball.jl")
# include("../../src/lcp.jl")
# include("../../src/solver.jl")

sys  = BouncingBall(Float64)
lcp  = ContactLCP.Lcp(Float64, sys)

X, t, Λn, Λt = ContactLCP.fulltimestep(lcp, sys.x0, Float64[]; totalTimeStep=1500)
plots(X, t, Λn, Λt)