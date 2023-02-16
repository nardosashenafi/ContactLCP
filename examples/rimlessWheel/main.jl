using ContactLCP
using LaTeXStrings
using PyPlot
using MeshCat
using GeometryBasics
using CoordinateTransformations
using ColorTypes
using Blink
using Rotations
using DiffEqFlux

include("dynamics.jl")
include("../../src/lcp.jl")
include("../../src/solver.jl")

# window = Window()
# vis = Visualizer()
# open(vis, window)

param = Float32[30.0, 5.0]
# param = Float32[0.0, 0.0]

sys  = RimlessWheel()
x0   = Float32.(initialState(pi, -0.5f0, 0.0f0, 0.0f0))
# x0 = Float32(initialState(sys, rand(-sys.α+0.1:0.05:0.0), rand(-2.0:0.05:-0.5), 0.0f0, 0.0f0))
# x0 = Float32([0.0, 0.26, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0])
lcp  = ContactLCP.Lcp(Float32, sys)

X, t, Λn, Λt = ContactLCP.fulltimestep(lcp, x0, param; expert=true, Δt = 0.001f0, totalTimeStep = 5000)
plots(X, t, Λn, Λt)

