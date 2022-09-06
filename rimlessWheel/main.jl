using LaTeXStrings, PyPlot
using JuMP, PATHSolver, Mosek, MosekTools, DiffOpt
using LinearAlgebra
import ChainRulesCore, ForwardDiff

include("dynamics.jl")

Δt = 0.001; totalTimeStep = 1000
θ0 = Float64[0.2]


sys  = RimlessWheel(Float64, Δt, totalTimeStep)
cm  = ContactMap(Float64, sys, θ0)
x0 = [0.2, lcp.sys.γ, -2.0, 0.0]

