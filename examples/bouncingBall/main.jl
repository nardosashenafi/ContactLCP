
using ContactLCP
using LaTeXStrings
using PyPlot

include("bouncing_ball.jl")


sys  = BouncingBall(Float32)
lcp  = ContactLCP.Lcp(Float32, sys)

X, t, Λn, Λt = ContactLCP.fulltimestep(lcp, sys.x0, Float32[]; totalTimeStep=1500)
plots(X, t, Λn, Λt)