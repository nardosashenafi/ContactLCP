
using ContactLCP
using LaTeXStrings
ENV["MPLBACKEND"]="tkagg"
using PyPlot

include("bouncing_ball.jl")


sys  = BouncingBall(Float64)
lcp  = ContactLCP.Lcp(sys, Float64)

X, t, Λn = ContactLCP.fulltimestep(lcp, Float64)
plots(X, t, Λn)
