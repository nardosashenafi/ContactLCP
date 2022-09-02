using ContactLCP
using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra

include("dynamics.jl")

Δt = 0.001; totalTimeStep = 500
θ0 = Float64[0.2]
x0 = [pi/2, 0.2, 0.0, 0.0 ]

sys  = RimlessWheel(Float64, Δt, totalTimeStep)
lcp  = ContactLCP.Lcp(Float64, sys, θ0)


function ContactLCP.fulltimestep(lcp::ContactLCP.Lcp, x0::Vector{T}, θ::Vector{T}) where {T<:Real}

    X       = Vector{Vector{Float64}}()
    Λn      = Vector{Vector{Float64}}()
    t       = Vector{T}()

    if isempty(x0)
        x = deepcopy(lcp.sys.x0)
    else
        x = deepcopy(x0)
    end

    X       = push!(X, x)
    t       = push!(t, T(0.0))

    for i in 1:lcp.sys.totalTimeStep
        if x[2] <= -lcp.sys.α
            x[2] = lcp.sys.α
        end

        x, λn = ContactLCP.oneTimeStep(lcp, x, θ)
        push!(X, x)
        push!(Λn, λn)
        push!(t, t[end]+lcp.sys.Δt)
    end

    return X, t, Λn
end


X, t, Λn = ContactLCP.fulltimestep(lcp, x0, θ0)
plots(X, t, Λn)