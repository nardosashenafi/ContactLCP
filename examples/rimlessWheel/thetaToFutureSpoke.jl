using ContactLCP
using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra

include("dynamics.jl")

Δt = 0.001; totalTimeStep = 1500
θ0 = Float64[0.2]


sys  = RimlessWheel(Float64, Δt, totalTimeStep)
lcp  = ContactLCP.Lcp(Float64, sys, θ0)
x0 = [lcp.sys.γ, 0.2, 0.0, -2.0 ]

function wrapPendulum(x)
    x[1] = atan(tan(x[1]))
    return x
end

function wrapWheel(lcp::ContactLCP.Lcp, x, λn )

    # if gap(lcp.sys, x) 
    if λn[1] > 0.0
        x[2] = -x[2]
    end

    return x
end

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
       
        x, λn = ContactLCP.oneTimeStep(lcp, x, θ)
        x = wrapWheel(lcp, x, λn)
        x = wrapPendulum(x)
        push!(X, x)
        push!(Λn, λn)
        push!(t, t[end]+lcp.sys.Δt)
    end

    return X, t, Λn
end

X, t, Λn = ContactLCP.fulltimestep(lcp, x0, θ0)
plots(X, t, Λn)