using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra
using ForwardDiff
using Flux

include("../dynamics.jl")
include("../../../src/lcp.jl")
include("../../../src/solver.jl")

sys  = RimlessWheel()
lcp  = Lcp(Float32, sys)
Δt   = 0.0005f0

function stateAndForcesWithNoise(lcp::Lcp, x, sysParam, controlParam::AbstractArray{T}; Δt = 0.001f0, kwargs...) where {T<:Real}

    qA, uA  = lcp.sys(x)
    qM      = qA + 0.5f0*Δt*uA

    x_mid   = vcat(qM, uA)
    gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold = sysAttributes(lcp, x_mid, sysParam, controlParam; kwargs...)
    contactIndex, _ = checkContact(x_mid, gn, gThreshold, k)
    λn, λt, λR  = solveLcp(gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold, x_mid; Δt=Δt)

    ######add noise on the angle of contact force to simulate uneven/uncertain terrain
    Wn_noise        = deepcopy(Wn)
    σx_Wn_noise     = 0.3f0 
    μx_Wn_noise     = 0.5f0 .* Wn_noise[1, contactIndex]   
    noise_nx        = σx_Wn_noise .* randn(length(contactIndex)) .+ μx_Wn_noise
    Wn_noise[1, contactIndex] += noise_nx 

    σy_Wn_noise  = 0.3f0 .* Wn_noise[2, contactIndex]
    μy_Wn_noise  = 0.5f0 .* Wn_noise[2, contactIndex]  
    noise_ny     = abs.(σy_Wn_noise .* randn(length(contactIndex)) .+ μy_Wn_noise)      #don't change the sign of the y component of the contact force in order to prevent it from penetrating into the ground
    Wn_noise[2, contactIndex] += noise_ny 

    ####complete integration
    uE = M\((Wn_noise - Wt*diagm(0 => μ))*λn + Wt*λR + h*Δt) + uA
    qE = qM + 0.5f0*Δt*uE

    return vcat(qE,uE), λn, λt
end

function oneTimeStepWithNoise(lcp::Lcp, x, sysParam, controlParam::AbstractArray{T}; Δt = 0.001f0, kwargs...) where {T<:Real}
    x2, _, _ = stateAndForcesWithNoise(lcp, x, sysParam, controlParam; Δt = Δt, kwargs...)
    return x2
end

oneStep(x, sysParam, controlParam; kwargs...) = oneTimeStepWithNoise(lcp, x, sysParam, controlParam; kwargs...)

function trajectory(x0, sysParam, controlParam::Vector{T}; expert=false, Δt = Δt, totalTimeStep = 1000) where {T<:Real}

    X       = Vector{Vector{T}}(undef, totalTimeStep)
    x       = deepcopy(x0)
 
    for i in 1:totalTimeStep
        x       = oneStep(x, sysParam, controlParam; Δt=Δt, expert=expert)
        X[i]    = x
    end
    # isStumbling(X) ? println("We got stumbling! Beware of gradient jumps") : nothing

    return X
end

function sampleSystemParam()
    γi = rand(0.0f0:0.01:30.0f0*pi/180f0)
    return γi
end

function sampleInitialStates(sysParam, controlParam::Vector{T}, sampleNum; α=α, totalTime=1000) where {T<:Real}

    sampleTrajectories = Vector{Vector{Vector{T}}}()

    for i in 1:2
        w  = rand(getq(controlParam))
        x0 = Float32.(initialState(rand(pi-α:0.05:pi+α), 
                                    rand(-3.0:0.05:-0.5), 
                                    0.0, 
                                    rand(-1.0:0.1:1.0)))


        S = trajectory(x0, sysParam, w; totalTimeStep=totalTime)
        push!(sampleTrajectories, S)
    end

    #sample 10 initial states from the long trajectories
    X0 = Vector{Vector{T}}()

    for i in 1:sampleNum
        if rand() < 0.5 
            push!(X0, rand(rand(sampleTrajectories)))
        else
            x0 = Float32.(initialState(rand(pi-α:0.05:pi+α), 
                                        rand(-3.0:0.05:-0.5), 
                                        0.0, 
                                        rand(-1.0:0.1:1.0)) )
            push!(X0, x0)
        end 
    end

    return X0

end