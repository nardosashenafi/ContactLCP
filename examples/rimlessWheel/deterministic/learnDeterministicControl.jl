
using DiffEqFlux
using MLBasedESC
using Statistics
using ProgressMeter, BSON
using MeshCat, GeometryBasics, CoordinateTransformations, ColorTypes, Blink, Rotations

include("deterministicTrainingHelpers.jl")

x0             = Float32.(initialState(pi, -1.0f0, 0.0f0, 0.0f0))
param_expert   = Float32[30.0, 5.0]
# unn            = FastChain(FastDense(6, 8, elu), 
#                             FastDense(8, 5, elu),
#                             FastDense(5, 1))
# ps             = 0.1f0*randn(Float32, DiffEqFlux.paramlength(unn))

Hd              = FastChain(FastDense(6, 8, elu), 
                  FastDense(8, 7, elu),
                  FastDense(7, 1))
const N         = 6
npbc            = MLBasedESC.NeuralPBC(N, Hd)
const satu      = 2.0f0

function extractStumbling(X)
    ind = findfirst(x -> x >= -0.01, getindex.(X, 8))
    if !isnothing(ind)
        return X[1:ind] 
    else
        return X
    end
end

function trajLoss(X0, param::Vector{T}; totalTime=500) where {T<:Real}

    l = 0.0f0
    for x0 in X0
        S       = trajectory(x0, param; totalTimeStep=totalTime)
        # S1      = extractStumbling(S)
        l       += hipSpeedLoss(S)
    end

    return sum(l)/length(X0)
end

function testControl(X0, ps, grad, fig1; timeSteps=5000)

    S       = trajectory(X0[1], ps; totalTimeStep = timeSteps);
    loss    = trajLoss(X0, ps; totalTime=timeSteps)

    plots(S, ps, fig1)

    println("loss = ", round(loss, digits=4) , " | grad = ", mean(grad) )
    # println("loss = ", round(l1(param), digits=4),  " | p = ", round.(param, digits=4), " | hip speed = ", round.(mean(-cm.sys.l1 .* cos.(getindex.(Z, 1)) .* getindex.(Z, 3)), digits=4) )
end

############################################################
#control to hip speed

function hipSpeedLoss(Z; gThreshold=gThreshold, k=k, α=α)

    #loss of one trajectory
    xd_dot  = 1.0f0

    loss = 0.0f0
    ki = 0

    for z in Z
        gn, _, _ = gap(z)
        contactIndex, _ = checkContact(z, gn, gThreshold, k)

        #loss between desired ẋ and actual ẋ should not be computed using getindex.(Z, 1) because this state does not directly depend on the control.
        #Using getindex.(Z, 1) as ẋ in auto-diff gives zero gradient
        if !isempty(contactIndex)
            ki = contactIndex[1] - 1
        else
            ki = spokeNearGround(gn) - 1
        end
        error = xd_dot - (l1 * cos(z[4] + 2*α*ki) * z[8])
        loss += 30.0f0*dot(error, error) + 3.0f0*dot(z[7], z[7])
    end

    return 1.0f0/length(Z)*loss
end

function controlToHipSpeed(;T=Float32)
    
    opt                 = Adam(0.005)

    counter             = 0
    X0                  = Vector{Vector{T}}()
    fig1                = plt.figure()
    param               = 0.3f0*randn(Float32, DiffEqFlux.paramlength(Hd)+N)
    param[end-N+1:end]  = 0.1f0*rand(N)
    minibatchSize       = 2
                            
    test(X, θ, grad)  = testControl(X, θ, grad, fig1; timeSteps=8000)

    @showprogress for i in 1:5000
        X0    = sampleInitialStates(param, minibatchSize; totalTime=7000)
        println("X0 = ", X0)

        l(θ)  = trajLoss(X0, θ; totalTime=4000)
        lg    = ForwardDiff.gradient(l, param)

        if counter > 4
            println(" ")
            test([X0[1]], param, lg)
            BSON.@save "./saved_weights/deterministic_hardware_6-8-8-7-7-1_elu.bson" param
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
#################################################################

#tracking expert trajectory
################################################################

function lossToExpert(lcp, param, x0; totalTime = 500)

    S, λ, _ = rwTrajectory(lcp, x0, param; totalTimeStep=totalTime)
    θ, θdot, ϕ, ϕdot = getindex.(S, 4), getindex.(S, 8), getindex.(S, 3), getindex.(S, 7)
    Sd, λd, _ = rwTrajectory(lcp, x0, param_expert; totalTimeStep=totalTime, expert=true); 
    θd, θdotd, ϕd, ϕdotd = getindex.(Sd, 4), getindex.(Sd, 8), getindex.(Sd, 3), getindex.(Sd, 7)

    return 500.0f0/length(S) * ( 
            2.0f0*dot(θd - θ , θd - θ) +
            2.0f0*dot(θdotd - θdot , θdotd - θdot) + 
            0.0f0*dot(ϕd - ϕ , ϕd - ϕ) + 
            0.0f0*dot(ϕdotd - ϕdot , ϕdotd - ϕdot) +  
            0.0f0*dot(λd - λ, λd - λ) )
end

function trackExpert(lcp::Lcp, x0::Vector{T}) where {T<:Real}

    param       = Float32[25.0, 2.0]
    opt         = Adam(0.01)

    l1          = Inf
    counter     = 0
    X0          = Vector{Vector{T}}()
    minibatchSize = 4

    for i in 1:400
        while isempty(X0)
            X0 = sampleInitialStates(param, minibatchSize; totalTime=1000)
        end
        for xi in X0
            l(θ)  = lossToExpert(lcp, θ, xi; totalTime = 500)
            lg    = ForwardDiff.gradient(l, param)

            l1    = l(param)
            if counter > 20
                println("loss = ", l1, " | grad = ", lg, " | param = ", param)
                counter = 0
            end
            counter += 1

            Flux.update!(opt, param, lg)
        end
    end
    println("param_expert= ", param_expert, " param_learned = ", param )
    return param
end

function gradientTest1(lcp, z, ps)

    #test the gradient of xdot wrt torque inputs
    function takeGradOf(lcp, z, param)
        S, λ, _ = rwTrajectory(lcp, z, param; totalTimeStep=500)
        # l = hipSpeedLoss(lcp, S)
        # u = control(z, param)
        # return l
        return sum(getindex.(S, 4))
    end

    l(θ) = takeGradOf(lcp, z, θ)
    lg   = ForwardDiff.gradient(l, ps)

end