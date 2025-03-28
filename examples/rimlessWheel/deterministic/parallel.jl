
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
                  FastDense(8, 5, elu),
                  FastDense(5, 1))
const N         = 6
npbc            = MLBasedESC.NeuralPBC(N, Hd)
const satu      = 3.0f0

function trajLoss(X0, param::Vector{T}; totalTime=500) where {T<:Real}

    l = Threads.Atomic{T}()
    Threads.@threads for x0 in X0
        lmin = oneBatch(x0, param;totalTimeStep=totalTime)
        Threads.atomic_add!(l, lmin)
    end

    return l.value/length(X0)
end

function oneBatch(x0, param;totalTimeStep=1500, α=α)
    S  = trajectory(x0, param; totalTimeStep=totalTimeStep)
    l  = hipSpeedLoss(S)

    # #add cost on contact frequency
    # β       = 0.2f0
    # xd_dot  = 0.5f0
    # freq_d  = 1.0f0/(2.0f0*l1*sin(α)/xd_dot)

    # θ             = getindex.(S, 4)
    # xdot          = getindex.(S, 5)
    # xddot_avg     = 1.0f0/length(θ)*sum(xdot)
    # strikePeriod  = 2.0f0*l1*sin(α)/xddot_avg  #TODO: assumes  if the spoke in contact
    # f             = 1.0f0/strikePeriod
    # if f >= (1+β)*freq_d
    #     l += 0.1f0*(f - (1+β) .* freq_d)
    # end
    return l
end

function testControl(X0, ps, grad, fig1; timeSteps=5000)

    S       = trajectory(X0[1], ps; totalTimeStep = timeSteps);
    loss    = trajLoss(X0, ps; totalTime=timeSteps)

    plots(S, fig1)

    println("loss = ", round(loss, digits=4) , " | grad = ", mean(grad) )
    # println("loss = ", round(l1(param), digits=4),  " | p = ", round.(param, digits=4), " | hip speed = ", round.(mean(-cm.sys.l1 .* cos.(getindex.(Z, 1)) .* getindex.(Z, 3)), digits=4) )
end

############################################################
#control to hip speed

function hipSpeedLoss(Z; gThreshold=gThreshold, k=k, α=α)

    #loss of one trajectory
    xd_dot  = 0.5f0

    loss = 0.0f0
    ki_prev = 0

    for z in Z
        gn, _, _ = gap(z)
        contactIndex, _ = checkContact(z, gn, gThreshold, k)

        #loss between desired ẋ and actual ẋ should not be computed using getindex.(Z, 1) because this state does not directly depend on the control.
        #Using getindex.(Z, 1) as ẋ in auto-diff gives zero gradient
        if !isempty(contactIndex)
            ki      = contactIndex[1] - 1
            loss    += xd_dot - (l1 * cos(z[4] + 2*α*ki) * z[8])
            ki_prev = ki
        else
            loss    += (xd_dot - (l1 * cos(z[4] + 2*α*ki_prev) * z[8]))    #the closest estimate to the hip speed when the spokes are not in contact
        end
    end

    lmag = dot(loss, loss)

    ϕdot = getindex.(Z, 7)
    lmag += 1.0f0*dot(ϕdot, ϕdot)

    # #add cost on contact frequency
    # β       = 0.2f0
    # freq_d  = 1.0f0/(2.0f0*lcp.sys.l1*sin(lcp.sys.α)/xd_dot)

    # xdot          = getindex.(X, 5)
    # xddot_avg     = 1.0f0/length(θ)*sum(xdot)
    # strikePeriod  = 2.0f0*l1*sin(α)/xddot_avg  #TODO: assumes  if the spoke in contact
    # f             = 1.0f0/strikePeriod
    # if f >= (1+β)*freq_d
    #     lmag += 0.1f0*(f - (1+β) .* freq_d)
    # end

    return 1.0f0/length(Z)*lmag
end

function gradient!(grad, X0, param::AbstractArray{T}; totalTimeStep = 1000) where {T<:Real}
    Threads.@threads for i in eachindex(X0)
        grad[i] = ForwardDiff.gradient((θ) -> oneBatch(X0[i], θ; totalTimeStep = totalTimeStep), param)
    end
end

function controlToHipSpeed(;T=Float32)
    
    opt                 = Adam(0.001)

    counter             = 0
    X0                  = Vector{Vector{T}}()
    fig1                = plt.figure()
    param               = 0.3f0*randn(Float32, DiffEqFlux.paramlength(Hd)+N)
    param[end-N+1:end]  = 0.1f0*rand(N)
    minibatchSize       = 2
    lg1                 = Vector{Vector{eltype(param)}}(undef, minibatchSize)
          
    test(X, θ, grad)  = testControl(X, θ, grad, fig1; timeSteps=10000)

    @showprogress for i in 1:5000
        println("Sampling")
        X0    = sampleInitialStates(param, minibatchSize; totalTime=3000)
        println("X0 = ", X0)
        # l(θ)  = trajLoss(X0, θ; totalTime=4000)
        # lg    = ForwardDiff.gradient(l, param)
        gradient!(lg1, X0, param; totalTimeStep = 2000) 
        lg   = mean(lg1)
        println("gradient complete")
        if counter > 4
            println(" ")
            test([X0[1]], param, lg)
            BSON.@save "./saved_weights/deterministic_hardware_6-8-8-5-5-1_elu.bson" param
            counter = 0
        end
        if any(isnan.(lg))
            println("NaN occurred")
            break
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