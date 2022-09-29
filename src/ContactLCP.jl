
module ContactLCP

using JuMP, PATHSolver
using LinearAlgebra
import ForwardDiff

include("lcp.jl")
include("solver.jl")

end # module