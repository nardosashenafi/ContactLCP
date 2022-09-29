using ContactLCP
using LaTeXStrings
using PyPlot

include("dynamics.jl")

Δt              = 0.001 
totalTimeStep   = 500

sys = RimlessWheel(Float64, Δt, totalTimeStep)
lcp  = Lcp(sys, Float64)

X, t, Λn, Λt = fulltimestep(lcp, sys.x0; totalTimeStep=1500)
plots(X, t)