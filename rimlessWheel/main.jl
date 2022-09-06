using LaTeXStrings, PyPlot
using JuMP, PATHSolver, Mosek, MosekTools, DiffOpt
using LinearAlgebra
import ChainRulesCore, ForwardDiff
using PyCall
@pyimport matplotlib.animation as anim

include("dynamics.jl")
include("contactMap.jl")

Δt = 0.001; totalTimeStep = 20000
θ0 = Float64[0.2]


sys  = RimlessWheel(Float64, Δt, totalTimeStep)
cm  = ContactMap(Float64, sys, θ0)
x0 = [0.2, cm.sys.γ, -2.0, 0.0]

X, t = fulltimestep(cm, x0, Float64[])
plots(X, t)

