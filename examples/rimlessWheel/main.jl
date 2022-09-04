using ContactLCP
ENV["MPLBACKEND"]="tkagg"
using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra

include("dynamics.jl")

Δt = 0.001; totalTimeStep = 1500
θ0 = Float64[0.2]


sys  = RimlessWheel(Float64, Δt, totalTimeStep)
lcp  = ContactLCP.Lcp(Float64, sys, θ0)
# x0 = [lcp.sys.γ, 0.2, 0.0, -2.0 ]
x0 = [0.5, 0.2, 0.0, -2.0 ]

function wrapPendulum(x)
    x[1] = atan(tan(x[1]))
    return x
end

function wrapWheel(lcp::ContactLCP.Lcp, x, λn)

    if (x[2] >= lcp.sys.α )   #moving backwards
        println("switching forward")
        x[2] = -lcp.sys.α
    elseif (x[2] <= -lcp.sys.α )   #moving forward
        println("switching backward")
        x[2] = lcp.sys.α
    end

    # g = gap(lcp.sys, x)
    # if λn[1] > 0.0
    #     x[2] = -x[2]
    # end

    return x
end

function impactMap(sys, ϕ)


    det = sys.I1*sys.I2 + sys.I1*sys.m2*sys.l2^2 + sys.I2*sys.mt*sys.l1^2 + 
            sys.m2*sys.l1^2*sys.l2^2*(sys.m1 + sys.m2*(sin(sys.α - ϕ))^2)
    ξ1 = 1/det * (sys.I1*sys.I2 + sys.I1*sys.m2*sys.l2^2) + 
        (sys.I2*sys.mt*sys.l1^2 + sys.m2*sys.l1^2*sys.l2^2*(sys.m1 +0.5* sys.m2))*cos(2sys.α) -
        1/2*sys.m2^2*sys.l1^2*sys.l2^2*cos(2ϕ)
    ξ2 = 1/det *(sys.m2*sys.l1*sys.l2*(sys.I1*(cos(sys.α - ϕ) - cos(sys.α + ϕ)) + 
        sys.mt*sys.l1^2*(cos(2sys.α)*cos(sys.α - ϕ) - cos(sys.α + ϕ))) )

    return [ξ1 0;
            ξ2 1]

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
       
        xpost, λn = ContactLCP.oneTimeStep(lcp, x, Float64[])
        if λn[1] > 0.0 
            #  xpost[2] = -xpost[2]
            #  xpost[3] = -xpost[3]
            #  xpost[4] = -xpost[4]

            # Complete algorithm hi-jacking
            xpost[2] = -xpost[2]
            xpost[3:4] = impactMap(sys, x[1])*x[3:4]
        end
        x = deepcopy(xpost)
        x = wrapPendulum(x)
        push!(X, x)
        push!(Λn, λn)
        push!(t, t[end]+lcp.sys.Δt)
    end

    return X, t, Λn
end

X, t, Λn = ContactLCP.fulltimestep(lcp, x0, θ0)
plots(X, t, Λn)
