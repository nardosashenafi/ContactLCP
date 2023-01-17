# using ContactLCP
using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra
using ForwardDiff

include("bouncing_ball.jl")
include("../../src/lcp.jl")
include("../../src/solver.jl")

Δt = 0.001f0; totalTimeStep = 1500
θ0 = Float32[0.2]
x0 = Float32[0.0,0.0002,0.2,-0.2]

sys  = BouncingBall(Float64)
lcp  = Lcp(Float64, sys)

function loss(lcp, θ, x0, Sθd, λθd; totalTimeStep = 1500)
    Sθ, λθ, _ = trajectory(lcp, x0, θ; totalTimeStep = totalTimeStep)
    return 50.0f0/length(Sθ) * ( 0.5f0*dot(Sθd - Sθ , Sθd - Sθ) + 0.5f0*dot(λθd - λθ, λθd - λθ) )
end

function oneStep(lcp::Lcp, x1, θm::AbstractArray{T}; Δt = 0.001f0) where {T<:Real}

    qA, uA  = lcp.sys(x1)
    qM      = qA + 0.5f0*Δt*uA

    x_mid   = vcat(qM, uA)
    # println("x = ", x_mid)
    gn, γn, γt, M, h, Wn, Wt = sysAttributes(lcp, x_mid, [])

    M = [θm[1] 0.0; 0.0 θm[1]]
    h = [0.0, -θm[1]*lcp.sys.g]
  
    λn, λt, λR  = solveLcp(gn, γn, γt, M, h, Wn, Wt, lcp.sys.ϵn, lcp.sys.ϵt, lcp.sys.μ, lcp.sys.gThreshold, x_mid; Δt=Δt) 
    # x2          = vcat(qM,uA)

    uE = M\((Wn - Wt*diagm(0 => lcp.sys.μ))*λn + Wt*λR + h*Δt) + uA
    qE = qM + 0.5f0*Δt*uE

    return vcat(qE,uE), λn, λt

end

function trajectory(lcp::Lcp, x0, θm::AbstractArray{T}; Δt = 0.001f0, totalTimeStep = 1500) where {T<:Real}

    X       = Vector{Vector{T}}(undef, totalTimeStep)
    Λn      = Vector{Vector{T}}(undef, totalTimeStep)
    Λt      = Vector{Vector{T}}(undef, totalTimeStep)
    x       = deepcopy(x0)

    for i in 1:totalTimeStep
        x, λn, λt  = oneStep(lcp, x, θm; Δt=Δt)
        X[i]    = x
        Λn[i]   = λn
        Λt[i]   = λt
    end

    return X, Λn, Λt
end

function solveM(lcp::Lcp, x0::Vector{T}) where {T<:Real}

    θm_actual   = [0.2]
    θm          = [0.9]
    η           = 0.2

    # x0 = [0.0, rand(range(0.2, 0.5, step=0.01)), 0.1, -0.1]
    Sθd, λθd, _ = trajectory(lcp, x0, θm_actual; totalTimeStep=100)
    l1 = Inf
    for i in 1:100
        l(θ)  = loss(lcp, θ, x0, Sθd, λθd; totalTimeStep=100)
        lg    = ReverseDiff.gradient(l, θm)
        θm    .-= η*lg
        l1    = l(θm)

        println("loss = ", l1, " | grad = ", lg, " | θm = ", θm)
    end
    println("θ_actual= ", θm_actual, " θ_learned = ", θm )
    return θm
end

# X, t, Λn, Λt = trajectory(lcp, x0, θ0)
# plots(X, t, Λn, Λt)
# solveM(lcp, x0)

function checkGradient(lcp::Lcp, x0::Vector{T}) where {T<:Real}

    θm_actual = [0.2]
    θ1        = [0.5]
    θ2        = θ1 .+ 0.001

    Sθd, _, λθd = trajectory(lcp, x0, θm_actual)
    l(θ)        = loss(lcp, θ, x0, Sθd, λθd)

    l1    = l(θ1)
    l2    = l(θ2)

    # lg      = ForwardDiff.gradient(l, θ1)
    l, back = Zygote.pullback(l, θ1)
    lg      = back(1)

    fd_l = (l2 - l1) / (θ2 - θ1)
    error_grad = abs.(fd_l - lg)

    println("error_grad = ", error_grad, " | fd_l = ", fd_l, " | lg = ", lg)
end
