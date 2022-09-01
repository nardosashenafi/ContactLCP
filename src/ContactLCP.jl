
module ContactLCP

using JuMP, PATHSolver, Mosek, MosekTools, DiffOpt
using LinearAlgebra
import ChainRulesCore, ForwardDiff

include("grads.jl")

end # module