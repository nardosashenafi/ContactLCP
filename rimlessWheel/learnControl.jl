
function sampleInitialStates(param; x0::Vector{T} = [0.2f0, 0.0f0, -2.0f0, 0.0f0]) where {T<:Real}

    sampleTrajectories = Vector{Vector{Float32}}()
    #generate 5 long trajectories
    while length(sampleTrajectories) < 2

        xi = Float32.([rand(-α:0.05:α), rand(-pi/2:0.1:pi/2), 
                    rand(-5.0:0.1:0.0), rand(-1.0:0.05:1.0)])

        X, tx = fulltimestep(xi, param; timeSteps=5000)
        sampleTrajectories = vcat(sampleTrajectories, X)

    end

    #sample 10 initial states from the long trajectories
    X0          = Vector{Vector{Float32}}()
    sampleNum   = 4
    for i in 1:sampleNum
        X0 = vcat(X0, [rand(rand(sampleTrajectories))])
    end

    if all(abs.(abs.(getindex.(X0, 1)) .- α) .< 0.001)
        for i in 0:Int(floor(sampleNum/2))-1
            X0[end-i] = rand(rand(sampleTrajectories)[1:3])
        end
    end
    return X0

end

function hipSpeedLoss(Z, tz)

    #loss of one trajectory
    β       = 0.3f0
    xd_dot  = 0.5f0
    freq_d  = 1.0f0/(2.0f0*l1*sin(α)/xd_dot)

    lnorm = 0.0f0

    for i in 1:length(Z)
        #each loop is one continuous part
        x       = Z[i]
        t       = tz[i]
        θ       = getindex.(x, 1)
        θdot    = getindex.(x, 3)
        l       = xd_dot .+ l1 .* cos.(θ) .* θdot 
        lnorm   += 1.0f0/length(θ)*dot(l, l) 

        #add cost on contact frequency
        if length(x) >= 2   #if a section of the trajectory contains only one state, strikePeriod will be 0 and frequency becomes inf
            xd_avg        = 1.0f0/length(θ)*sum(abs.(-l1 .* cos.(θ) .* θdot))
            strikePeriod  = 2.0f0*l1*sin(α)/xd_avg
            f             = 1.0f0/strikePeriod
            if f >= (1+β)*freq_d
                lnorm += 0.1f0*(f - (1+β) .* freq_d)
            end
        end
    end

    return 1.0f0/length(Z)*lnorm
end

function trajLoss(x0, param)

    # gn, M, B, C, G, τ, Ξ, contactIndex, Y, tY = Array.(allocateCuArrays(x0, param))
    # loss = fulltimestep!(x0, param, gn, M, B, C, G, τ, Ξ, contactIndex, Y, tY; timeSteps=1000)
    X, tx = fulltimestep(x0, param; timeSteps=1000)
    loss = hipSpeedLoss(X, tx)

    return loss
end

function gradient!(X0, param, grad)

    batchsize = length(X0)

    for i in 1:batchsize
        l(θ)   = trajLoss(CuArray(X0[i]), θ) 
        _, back = Zygote.pullback(l, param)
        grad[i] = back(1.0f0)
    end

end

function controlToHipSpeed(ps)

    counter             = 0
    opt                 = Flux.Adam(0.001)
    param               = deepcopy(ps)
    X0                  = Vector{Vector{Float32}}()
    fig1                = plt.figure()
 

    for i in 1:2000

        while isempty(X0)   #this check is necessary if stumbling is extracted. DAgger may return empty replay buffer
            X0 = sampleInitialStates(Array(ps))
        end

        grad = Vector{typeof(param)}(undef, length(X0))
        
        CUDA.@sync begin
            gradient!(X0, param, grad)
        end

        if counter > 1
            testControl(X0[end], param, grad, fig1; timeSteps=5000)
            counter = 0
        end

        Flux.update!(opt, param, mean(grad))
        counter += 1

        X0  = Vector{Vector{Float64}}()
    end

    return param
end