
using DiffEqFlux
using Statistics
include("trainingHelpers.jl")

Δt = 0.001f0; totalTimeStep = 1500
sys            = RimlessWheel(Float32)
x0             = initialState(sys, 0.0f0, -0.5f0, 0.0f0, 0.0f0)
param_expert   = Float32[30.0, 5.0]
lcp            = Lcp(Float32, sys)
unn            = FastChain(FastDense(6, 8, elu), FastDense(8, 1))
ps             = 0.05*randn(DiffEqFlux.paramlength(unn))
const satu     = 1.5

function trajLoss(lcp::Lcp, X0, param::Vector{T}; totalTime=500) where {T<:Real}

    l = 0.0
    for x0 in X0
        S, λ, _ = rwTrajectory(lcp, x0, param; totalTimeStep=totalTime)
        # X, tx   = extractStumbling(X, tx)
        l       += hipSpeedLoss(lcp, S)
        l       += 0.05/length(S)*sum([abs.(unn(inputLayer(s), param)[1]) for s in S])
    end

    return sum(l)/length(X0)
end


function testControl(lcp::Lcp, x0, ps, grad, fig1; timeSteps=1000)

    S, λ, _ = rwTrajectory(lcp, x0, ps; totalTimeStep = timeSteps);
    # loss    = hipSpeedLoss(lcp, S)
    # loss  += 1.0/length(S)*sum([abs.(unn(inputLayer(s), param)[1]) for s in S])

    loss = trajLoss(lcp, [x0], ps; totalTime=timeSteps)

    plots(lcp.sys, S, fig1)

    println("loss = ", round(loss, digits=4) , " | grad ", norm(grad, 2), " | hip speed = ", round.(mean(getindex.(S, 5)), digits=4) )
    # println("loss = ", round(l1(param), digits=4),  " | p = ", round.(param, digits=4), " | hip speed = ", round.(mean(-cm.sys.l1 .* cos.(getindex.(Z, 1)) .* getindex.(Z, 3)), digits=4) )
end

############################################################
#control to hip speed

function hipSpeedLoss(lcp::Lcp, X)

    #loss of one trajectory
    # β       = 0.2f0
    xd_dot  = 0.5f0
    # freq_d  = 1.0f0/(2.0f0*lcp.sys.l1*sin(lcp.sys.α)/xd_dot)

    θ       = getindex.(X, 4)
    θdot    = getindex.(X, 8)
    l       = xd_dot .+ lcp.sys.l1 .* cos.(θ) .* θdot 
    lnorm   = dot(l, l)

    #add cost on contact frequency
    # if length(x) >= 2   #if a section of the trajectory contains only one state, strikePeriod will be 0 and frequency becomes inf
    #     xd_avg        = 1/length(θ)*sum(x)
    #     strikePeriod  = 2*cm.sys.l1*sin(cm.sys.α)/xd_avg  #TODO: assumes  if the spoke in contact
    #     f             = 1/strikePeriod
    #     if f >= (1+β)*freq_d
    #         lnorm += 0.1f0*(f - (1+β) .* freq_d)
    #     end
    # end

    return 10.0f0/length(X)*lnorm
end

function controlToHipSpeed(lcp::Lcp, x0::Vector{T}, ps) where {T<:Real}
    
    opt         = Adam(0.001)

    l1          = Inf
    counter     = 0
    X0          = Vector{Vector{T}}()
    fig1        = plt.figure()
    param       = deepcopy(ps)

    for i in 1:1000
        while isempty(X0)
            X0 = sampleInitialStates(lcp, param; totalTime=1000)
        end

        l(θ)  = trajLoss(lcp, X0, θ;  totalTime=500)
        lg    = ForwardDiff.gradient(l, param)

        if counter > 10
            testControl(lcp, x0, param, lg, fig1)
            counter = 0
        end

        counter += 1

        Flux.update!(opt, param, lg)

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

    for i in 1:400
        while isempty(X0)
            X0 = sampleInitialStates(lcp, param; totalTime=1000)
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