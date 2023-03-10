using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra
using ForwardDiff
using Flux

include("../dynamics.jl")
include("../../../src/lcp.jl")
include("../../../src/solver.jl")

sys  = RimlessWheel()
lcp  = Lcp(Float32, sys)
Δt   = 0.0002f0


function unstackParams(param)
    μ_param = @view param[1:paramNum]
    σ_param = @view param[paramNum+1:end]

    return μ_param, σ_param
end

function stateAndForcesWithNoise(lcp::Lcp, x, sysParam, controlParam::AbstractArray{T}; Δt = 0.001f0, kwargs...) where {T<:Real}

    qA, uA  = lcp.sys(x)
    qM      = qA + 0.5f0*Δt*uA

    x_mid   = vcat(qM, uA)
    gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold = sysAttributes(lcp, x_mid, sysParam, controlParam; kwargs...)
    λn, λt, λR  = solveLcp(gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold, x_mid; Δt=Δt)

    ####complete integration
    uE = M\((Wn - Wt*diagm(0 => μ))*λn + Wt*λR + h*Δt) + uA
    qE = qM + 0.5f0*Δt*uE

    return vcat(qE,uE), λn, λt
end

function oneTimeStepWithNoise(lcp::Lcp, x, sysParam, controlParam::AbstractArray{T}; Δt = 0.001f0, kwargs...) where {T<:Real}
    x2, _, _ = stateAndForcesWithNoise(lcp, x, sysParam, controlParam; Δt = Δt, kwargs...)
    return x2
end

function oneTimeStepWithNoise(lcp::Lcp, x, sysParam, u::T; Δt = 0.001f0, kwargs...) where {T<:Real}

    qA, uA  = lcp.sys(x)
    qM      = qA + 0.5f0*Δt*uA

    x_mid   = vcat(qM, uA)
    gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold = sysAttributes(lcp, x_mid, sysParam, u; kwargs...)
    λn, λt, λR  = solveLcp(gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold, x_mid; Δt=Δt)

    ####complete integration
    uE = M\((Wn - Wt*diagm(0 => μ))*λn + Wt*λR + h*Δt) + uA
    qE = qM + 0.5f0*Δt*uE

    return vcat(qE,uE)
end

oneStep(x, sysParam, controlParam::AbstractArray{T}; kwargs...) where {T<:Real} = oneTimeStepWithNoise(lcp, x, sysParam, controlParam; kwargs...)
oneStep(x, sysParam, u::T; kwargs...) where {T<:Real} = oneTimeStepWithNoise(lcp, x, sysParam, u; kwargs...)

function trajectory(x0, r, controlParam::Vector{T}; expert=false, Δt = Δt, totalTimeStep = 1000) where {T<:Real}

    X           = Vector{Vector{T}}(undef, totalTimeStep)
    obstacles   = Vector{Vector{Vector{T}}}(undef, totalTimeStep)
    x           = deepcopy(x0)
    sysParam    = createUnevenTerrain(x, r)
    θi          = x[4]

    for i in 1:totalTimeStep

        if abs(x[4] - θi) > 2pi - 2α
            sysParam = createUnevenTerrain(x, r)
            θi = x[4]
        end
        x       = oneStep(x, sysParam, controlParam; Δt=Δt, expert=expert)
        X[i]    = x
        obstacles[i] = sysParam
    end

    return X, obstacles
end

function hipSpeedLoss(Z, obstacles; gThreshold=gThreshold, k=k, α=α)

    #loss of one trajectory
    xd_dot  = 0.5f0
    loss = 0.0f0
    ki = 0

    for i in eachindex(Z)
        z = Z[i]
        gn, _, _ = gap(z, obstacles[i])
        contactIndex, _ = checkContact(z, gn, gThreshold, k)

        #loss between desired ẋ and actual ẋ should not be computed using getindex.(Z, 1) because this state does not directly depend on the control.
        #Using getindex.(Z, 1) as ẋ in auto-diff gives zero gradient
        if !isempty(contactIndex)
            ki = contactIndex[1] - 1
        else
            ki = spokeNearGround(gn) - 1
        end
        error = xd_dot - (l1 * cos(z[4] + 2*α*ki) * z[8])
        loss += 10.0f0*dot(error, error) + 0.1f0*dot(z[7], z[7])
    end

    return 1.0f0/length(Z)*loss
end

function isStumbling(x)
    x[8] >= -0.01
end

function sampleInitialStates(controlParam::Vector{T}, sampleNum; α=α, totalTime=1000) where {T<:Real}

    sampleTrajectories = Vector{Vector{Vector{T}}}()
    rmax = 0.03f0

    w     = rand(getq(controlParam))
    x0, r = initialStateWithBumps(
                rand(pi-α+0.2f0:0.05f0:pi+α-0.2f0), 
                rand(-3.0f0:0.05f0:0.0f0), 
                0.0f0, 
                rand(-1.0f0:0.1f0:1.0f0), rmax)


    S, _ = trajectory(x0, r, w; totalTimeStep=totalTime)
    push!(sampleTrajectories, S)

    #sample 10 initial states from the long trajectories
    X0 = Vector{Vector{T}}()
    R = Vector{Vector{T}}()

    for i in 1:sampleNum
        if rand() < 0.5 
            push!(X0, rand(rand(sampleTrajectories)))
            push!(R, r)
        else
            x0, r1 = initialStateWithBumps(
                    rand(pi-α+0.2f0:0.05f0:pi+α-0.2f0), 
                    rand(-3.0f0:0.05f0:0.0f0), 
                    0.0f0, 
                    rand(-1.0f0:0.1f0:1.0f0), rmax)

            push!(X0, x0)
            push!(R, r1)
        end 
    end
    #extract stumbling
    for i in eachindex(X0)
        if isStumbling(X0[i])
            X0[i], R[i] = initialStateWithBumps(
                rand(pi-α+0.2f0:0.05f0:pi+α-0.2f0), 
                rand(-3.0f0:0.05f0:0.0f0), 
                rand(-pi/4.0f0:0.05f0:pi/4.0f0), 
                rand(-1.0f0:0.1f0:1.0f0), rmax)
        end

    end
    return X0, R

end

function hasconverged(x0, r, param, elbo_data, i)

    filtered_elbo = filter_elbo(elbo_data[end-4:end])

    println("elbo = ", round(elbo_data[end], digits=4))

    w                   = rand(getq(param))
    X, obstacles        = trajectory(x0, r, w; totalTimeStep=5000)
    loss                = hipSpeedLoss(X, obstacles)
    plots(X, fig1)
    ax2.plot(i, filtered_elbo[end], marker=".", color="k") 
    BSON.@save "./saved_weights/RW_bayesian_6-8-8-5-5-1_elu.bson" param
    println("loss = ", round(loss, digits=4) , " | hip speed = ", round.(mean(getindex.(X, 5)), digits=4) )

    return false
end

function filter_elbo(Δelbo_vec)
    DSP.filt(DSP.digitalfilter(responsetype, designmethod), Δelbo_vec)
end

function plots(Z, fig1)
    PyPlot.figure(1)
    fig1.clf()
    subplot(2, 2, 1)
    plot(getindex.(Z, 3), getindex.(Z, 7))
    scatter(Z[end][3], Z[end][7])
    ylabel(L"\dot{\phi} [rad/s]", fontsize=15)
    subplot(2, 2, 2)
    # θ, impactIndex = spokeInContact(Z)
    plot(getindex.(Z, 4), getindex.(Z, 8))
    scatter(Z[end][4], Z[end][8])
    ylabel(L"\dot{\theta} [rad/s]", fontsize=15)
    subplot(2, 2, 3)
    plot(getindex.(Z, 5))
    ylabel("vx [m/s]", fontsize=15)
    println("Average hip speed = ", mean(getindex.(Z, 5)))
end

function integrateMarginalization(x0, r, controlParam::Vector{T}, sampleNum; expert=false, Δt = Δt, totalTimeStep = 1000) where {T<:Real}

    X           = Vector{Vector{T}}(undef, totalTimeStep)
    obstacles   = Vector{Vector{Vector{T}}}(undef, totalTimeStep)
    sysParam    = createUnevenTerrain(x0, r)
    x           = deepcopy(x0)
    θi          = x[4]

    for i in 1:totalTimeStep
        if abs( x[4] - θi) > 2pi - 2α
            sysParam = createUnevenTerrain(x, r)
            θi =  x[4]
        end
        u       = marginalize(x, controlParam; sampleNum=sampleNum)
        x       = oneStep(x, sysParam, u; Δt=Δt, expert=expert)
        X[i]    = x
        obstacles[i] = sysParam
    end

    return X, obstacles
    
end

function control(z, u::T; expert=false) where {T<:Real}
    q, v = parseStates(z)

    if expert #working expert controller
        @assert length(θp) == 2
        return -θp[1]*(q[3]-0.65f0) - θp[2]*v[3]
    else    
        return clamp(u, -satu, satu)
    end
end

function genForces(z, u::T; expert=false) where {T<:Real}

    q, v = parseStates(z)
    x, y, ϕ, θ = q
    xdot, ydot, ϕdot, θdot = v
    #h = Bu - C qdot - G
    B = [0.0f0, 0.0f0, 1.0f0, -1.0f0]

    C = [0.0f0 0.0f0 -m2*l2*sin(ϕ)*ϕdot 0.0f0;
        0.0f0 0.0f0 m2*l2*cos(ϕ)*ϕdot 0.0f0;
        -m2*l2*sin(ϕ)*ϕdot m2*l2*cos(ϕ)*ϕdot 0.0f0 0.0f0;
        0.0f0 0.0f0 0.0f0 0.0f0]

    G = [-mt*g*sin(γ), 
        mt*g*cos(γ), 
        m2*g*l2*sin(ϕ - γ),
        0.0f0]

    return B*control(z, u; expert=expert) - C*v - G
end

function (sys::RimlessWheel)(z, sysParam, u::T; ϵn=ϵn_const, ϵt=ϵt_const, μ=μ_const, gThreshold=gThreshold, expert=false) where {T<:Real}
    q, v = parseStates(z)
    x, y, ϕ, θ = q
    gn, Wn, Wt = gap(z, sysParam)
    γn  = vnormal(Wn, v)
    γt  = vtang(Wt, v)
    M   = massMatrix(ϕ)
    h   = genForces(z, u; expert=expert)

    return gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold
end
