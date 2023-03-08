using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra
using ForwardDiff
using Flux
using DiffEqFlux
using MLBasedESC
using Statistics
using ProgressMeter, BSON
using MeshCat, GeometryBasics, CoordinateTransformations, ColorTypes, Blink, Rotations


include("../dynamics.jl")
include("../../../src/lcp.jl")
include("../../../src/solver.jl")

sys  = RimlessWheel()
lcp  = Lcp(Float32, sys)
Δt   = 0.0002f0

x0             = Float32.(initialState(pi, -1.0f0, 0.0f0, 0.0f0))
param_expert   = Float32[30.0, 5.0]

Hd              = FastChain(FastDense(6, 6, elu), 
                  FastDense(6, 4, elu),
                  FastDense(4, 1))
const N         = 6
npbc            = MLBasedESC.NeuralPBC(N, Hd)
const satu      = 1.5f0


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

oneStep(x, sysParam, controlParam; kwargs...) = oneTimeStepWithNoise(lcp, x, sysParam, controlParam; kwargs...)

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

function trajLoss(X0, R, param::Vector{T}; totalTime=500) where {T<:Real}

    l11 = 0.0f0
    for i in eachindex(X0)
        S, obstacles = trajectory(X0[i], R[i], param; totalTimeStep=totalTime)
        l11 += hipSpeedLoss(S, obstacles)
    end

    return sum(l11)/length(X0)
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

function sampleInitialStatesUnevenTerrain(controlParam::Vector{T}, r_const, sampleNum; α=α, totalTime=1000) where {T<:Real}

    sampleTrajectories = Vector{Vector{Vector{T}}}()

    x0, r = initialStateWithBumps(
                rand(pi-α+0.2f0:0.05f0:pi+α-0.2f0), 
                rand(-3.0f0:0.1f0:0.0f0), 
                0.0f0, 
                rand(-1.0f0:0.1f0:1.0f0), r_const)


    S, _ = trajectory(x0, r, controlParam; totalTimeStep=totalTime)
    push!(sampleTrajectories, S)

    #sample sampleNum initial states from the long trajectories
    X0 = Vector{Vector{T}}()
    R = Vector{Vector{T}}()

    for i in 1:sampleNum
        if rand() < 0.5 
            push!(X0, rand(rand(sampleTrajectories)))
            push!(R, r)
        else
            x0, r1 = initialStateWithBumps(
                    rand(pi-α+0.2f0:0.05f0:pi+α-0.2f0), 
                    rand(-3.0f0:0.1f0:0.0f0), 
                    0.0f0, 
                    rand(-1.0f0:0.1f0:1.0f0), r_const)

            push!(X0, x0)
            push!(R, r1)
        end 
    end
    #extract stumbling
    for i in eachindex(X0)
        if isStumbling(X0[i])
            X0[i], R[i] = initialStateWithBumps(
                            rand(pi-α+0.2f0:0.05f0:pi+α-0.2f0), 
                            rand(-3.0f0:0.1f0:0.0f0), 
                            rand(-pi/4.0f0:0.05f0:pi/4.0f0), 
                            rand(-1.0f0:0.1f0:1.0f0), r_const)
        end

    end
    return X0, R

end

function testControl(X0, R, ps, grad, fig1; timeSteps=5000)

    S,_       = trajectory(X0[1], R[1], ps; totalTimeStep = timeSteps);
    loss    = trajLoss(X0, R, ps; totalTime=timeSteps)

    plots(S, fig1)

    println("loss = ", round(loss, digits=4) , " | grad = ", maximum(grad) )
    # println("loss = ", round(l1(param), digits=4),  " | p = ", round.(param, digits=4), " | hip speed = ", round.(mean(-cm.sys.l1 .* cos.(getindex.(Z, 1)) .* getindex.(Z, 3)), digits=4) )
end

function controlToHipSpeedUnvevenTerrain(;T=Float32)
    
    opt                 = Adam(0.005)

    counter             = 0
    X0                  = Vector{Vector{T}}()
    fig1                = plt.figure()
    param               = 0.3f0*randn(Float32, DiffEqFlux.paramlength(Hd)+N)
    param[end-N+1:end]  = 0.1f0*rand(N)
    minibatchsize       = 2
                            
    test(X, R, θ, grad)  = testControl(X, R, θ, grad, fig1; timeSteps=8000)
    _, r_const           = initialStateWithBumps(
                                rand(pi-α+0.2f0:0.05f0:pi+α-0.2f0), 
                                rand(-4.0f0:0.05f0:-1.0f0), 
                                0.0f0, 
                                rand(-1.0f0:0.1f0:1.0f0), 0.03f0)

    @showprogress for i in 1:5000
        X0, R = sampleInitialStatesUnevenTerrain(param, r_const, minibatchsize; 
                                                    totalTime=6000)
        println("X0 = ", X0)
        println("R = ", R)

        l(θ)  = trajLoss(X0, R, θ; totalTime=6000)
        lg    = ForwardDiff.gradient(l, param)

        if counter > 4
            println(" ")
            test([X0[1]], [R[1]], param, lg)
            BSON.@save "./saved_weights/unevenTerrain_deterministic_hardware_6-6-6-4-4-1_elu.bson" param
            counter = 0
        end
        if any(isnan.(lg))
            println("NaN occurred")
            continue
        end
        Flux.update!(opt, param, lg)
        counter += 1
    end

    return param
end