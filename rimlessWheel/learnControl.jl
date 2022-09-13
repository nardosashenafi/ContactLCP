
function interpolate_τ(limitcycle_x, limitcycle_t)
    limitcycle_θ    =  getindex.(limitcycle_x, 1)
    τ               = Interpolations.extrapolate(Interpolations.interpolate((-limitcycle_θ,), limitcycle_t, Gridded(Linear())), Line())
    return τ
end

function interpolate_θdot(limitcycle_x, limitcycle_t)
    limitcycle_θdot =  getindex.(limitcycle_x, 3)
    θdot_star       = Interpolations.extrapolate(Interpolations.interpolate((limitcycle_t,), limitcycle_θdot, Gridded(Linear())), Line())
    return θdot_star
end

function isStumbling(x)
    # length(x) <= 50
    (any(getindex.(x, 3) .> 0.0) || length(x) <= 2)     
    # false
end

function extractStumbling(X::Vector{Vector{Vector{T}}}, tx) where {T<:Real}
    Z       = Vector{Vector{Vector{T}}}()
    tz      = Vector{Vector{T}}()
    for (x, t) in collect(zip(X, tx))
        if !isStumbling(x)
            push!(Z, x)
            push!(tz, t)
        end
    end

    return Z, tz
end

function perpLoss(x0, τ, θdotstar, param::Vector{T}) where {T<:Real}

    X, tx = fulltimestep(cm, x0, param; timeSteps=500)
    loss  = 0.0
    Z     = Vector{Vector{Vector{T}}}()

    Z, tz = extractStumbling(X, tx)

    for x in Z
        θ           = getindex.(x, 1)
        θdot        = getindex.(x, 3)
        ϕdot        = getindex.(x, 4)
        τi          = τ.(-θ)
        θdotstar_i  = θdotstar.(τi)
        loss        += 1.0/length(θ)*(sum(map((xs, s) -> dot(xs-s, xs-s), θdot, θdotstar_i)) +
                        0.0*dot(ϕdot, ϕdot))

    end

    return loss/length(Z)
end

function hipSpeedLoss(Z, tz)

    #loss of one trajectory
    β       = 0.3
    xd_dot  = 0.5
    freq_d  = 1.0/(2*cm.sys.l1*sin(cm.sys.α)/xd_dot)

    lnorm = 0.0
    for (x, t) in collect(zip(Z, tz))
        #each loop is one continuous part
        θ       = getindex.(x, 1)
        θdot    = getindex.(x, 3)
        ϕdot    = getindex.(x, 4)
        l       = xd_dot .+ cm.sys.l1 .* cos.(θ) .* θdot 
        lnorm    += 1/length(θ)*dot(l, l) 

        #add cost on contact frequency
        # if length(x) >= 2
        #     strikePeriod  = t[end] - t[1] 
        #     f       = 1/strikePeriod
        #     if f >= (1+β)*freq_d
        #         lnorm += 1.0*(f - (1+β) .* freq_d)
        #     end
        # end
    end

    return 1.0/length(Z)*lnorm
end

# function trajLoss(cm, x0::Vector{T}, param) where {T<:Real}
function trajLoss(cm, x0, param::Vector{T}) where {T<:Real}

    X, tx   = fulltimestep(cm, x0, param; timeSteps=1000);
    X, tx   = extractStumbling(X, tx)

    return hipSpeedLoss(X, tx)
end

# function sampleInitialStates(x0::Vector{T}, param) where {T<:Real}
function sampleInitialStates(x0, param::Vector{T}) where {T<:Real}

    Z = Vector{Vector{Vector{T}}}()

    while isempty(Z)
        X, tx = fulltimestep(cm, x0, param; timeSteps=5000)
        Z, tz = extractStumbling(X, tx) #TODO: maybe slowly adding stumbling
    end

    X0 = Vector{Vector{T}}()
    sampleNum = 50
    for i in 1:sampleNum
        push!(X0, rand(rand(Z)))
    end
    return X0

end

function controlToHipSpeed(cm::ContactMap, ps)

    counter     = 0
    opt         = Flux.Adam(0.001)
    # param       = 0.2*rand(Lux.parameterlength(unn)) 
    param       = deepcopy(ps)
    fig1        = plt.figure()
    loss        = Inf
    X0          = Vector{Vector{Float64}}()
    # opt_state   = Optimisers.setup(opt, param)

    for i in 1:2000

        while isempty(X0)
            x0 = [rand(-cm.sys.α:0.05:cm.sys.α), rand(-pi/2:0.1:pi/2), rand(-5.0:0.1:-2.0), 0.0]
            X0 = sampleInitialStates(x0, param)
        end

        for xi in X0
            l1(θ)   = trajLoss(cm, xi, θ)
            grad    = ForwardDiff.gradient(l1, param)
            # pf_f = AD.value_and_pushforward_function(AD.ForwardDiffBackend(), l1, θ)
            # grad = Lux.gradient(l1, param)        #only uses Zygote, which does not like mutating arrays

            if counter > 30
                testControl(cm, x0, param, grad, fig1; timeSteps=5000)
                counter = 0
            end

            Flux.update!(opt, param, grad)
            # opt_state, param = Optimisers.update(opt_state, param, grad)
            counter += 1
            # loss = l1(param)
        end
        X0  = Vector{Vector{Float64}}()
    end

    return param
end

function testControl(cm, x0, param, grad, fig1; timeSteps=5000)
    X, tx   = fulltimestep(cm, x0, param; timeSteps=timeSteps);
    loss    = hipSpeedLoss(X, tx)
    Z       = reduce(vcat, X)

    Ze, tze = extractStumbling(X, tx)   #extracting stumbles to compute hip speed
    Ze      = reduce(vcat, Ze)

    fig1.clf()
    subplot(2, 1, 1)
    plot(getindex.(Z, 1), getindex.(Z, 3))
    subplot(2, 1, 2)
    plot(getindex.(Z, 2))

    println("loss = ", round(loss, digits=4) , " | hip speed = ", round.(mean(-cm.sys.l1 .* cos.(getindex.(Ze, 1)) .* getindex.(Ze, 3)), digits=4) )
    # println("loss = ", round(l1(param), digits=4),  " | p = ", round.(param, digits=4), " | hip speed = ", round.(mean(-cm.sys.l1 .* cos.(getindex.(Z, 1)) .* getindex.(Z, 3)), digits=4) )
end

function controlToLimitCycle(cm::ContactMap) 

    S, t        = limitCycle(cm)
    τ           = interpolate_τ(S, t)
    θdotstar    = interpolate_θdot(S, t)
    counter     = 0
    param       = deepcopy(controlparam)
    opt         = Adam(0.005)
    fig1        = plt.figure()
    loss        = Inf
    X0          = Vector{Vector{T}}()

    for i in 1:2000
        while isempty(X0)
            # x0 = [rand(-cm.sys.α:0.05:cm.sys.α), rand(0:0.1:pi/2), rand(-5.0:0.1:0.0), 0.0]
            X0 = sampleInitialStates(x0, param)
        end

        for xi in X0
            l1(θ)   = perpLoss(xi, τ, θdotstar , θ)
            grad    = ForwardDiff.gradient(l1, param)
            if counter > 10
                println("loss = ", l1(param), " θ = ", param)
                println("grad = ", grad)
                X, t    = fulltimestep(cm, xi, param; timeSteps=5000)
                Z       = reduce(vcat, X)
                fig1.clf()
                subplot(2, 1, 1)
                plot(getindex.(Z, 1), getindex.(Z, 3))
                plot(getindex.(S, 1), getindex.(S, 3))
                subplot(2, 1, 2)
                plot(getindex.(Z, 2))
                counter = 0
            end

            Flux.Optimise.update!(opt, param, grad)
            counter += 1
            loss = l1(param)
        end
    end

end


function fd(x0, τ, θdotstar, Δ, param1)
    param2 = param1 + [Δ, 0.0]
    l1 = perpLoss(x0, τ, θdotstar, param1) 
    l2 = perpLoss(x0, τ, θdotstar, param2)

    param3 = param1 + [0.0, Δ]
    l3 = perpLoss(x0, τ, θdotstar, param3)
    return [(l2 - l1) ./ (param2[1] - param1[1]), (l3 - l1) ./ (param3[2] - param1[2])]
end

function matchExpert(cm::ContactMap, x0::Vector{T}) where {T<:Real}

    function matchloss(x0::Vector{T}, S, param) where {T<:Real}
        X, t = fulltimestep(cm, x0, param)
    
        Q = diagm(0 => [2.0, 1.0, 0.5, 0.5])
        l = 5.0/length(X) * sum(map((x, s) -> dot(x-s, Q*(x-s)), X, S))
        return l
    end

    S, t    = fulltimestep(cm, x0, [20.0, 5.0])

    param   = [50.0, 10.0]
    l1(θ)   = matchloss(x0, S, θ)
    opt     = Adam(0.1)
    counter = 0

    while l1(param) > 0.001
        grad = ForwardDiff.gradient(l1, param)
        Flux.Optimise.update!(opt, param, grad)

        if counter > 50
            println("loss = ", l1(param), " grad = ", grad, " θ = ", param)
            counter = 0
        end
        counter += 1
    end

    return param
end
