include("EMSwingup.jl")

function sampleInitialState(param::Vector{T}; totalTimeStep = totalTimeStep, minibatch=4) where {T<:Real}
    
    x0  = T.([rand(-TRACK_LENGTH/2.0:0.01:TRACK_LENGTH/2.0), rand(0.0:0.005:2pi), 
                  rand(-0.5:0.001:0.5), rand(-3.0:0.001:3.0)])

    X           = trajectory(x0, param)
    samples     = Vector{Vector{T}}(undef, minibatch)
    cropBuffer  = 0.3f0

    for i in 1:minibatch
        if rand() > 0.2
            select     = rand(X)
            #It is highly likely that the initial parameter cause the cart to go far from x1 = 0
            #Hence sampling from this trajectory tends to overwhelm the reply buffer with large x1 values.
            #Clamping x1 is bad because it tends to always work under the edge cases
            #so simply randomize x1 when it goes too far
            if  !(-TRACK_LENGTH/2.0 <= select[1] <= TRACK_LENGTH/2.0)
                select[1]  = rand(-TRACK_LENGTH/2.0:0.01:TRACK_LENGTH/2.0) 
            end
            samples[i] = select
        else
            samples[i] =  [rand(-0.1:0.001:0.1), rand(-0.5:0.001:0.5), 
                            rand(-0.5:0.001:0.5), rand(-1.0:0.001:1.0)]
        end
    end

    return samples
end

function trajectory(xi, param::AbstractArray{T}; totalTimeStep = totalTimeStep) where {T<:Real}

    X = Vector{Vector{T}}(undef, totalTimeStep)
    X[1] = xi

    for i in 1:totalTimeStep-1
        u       = controlArray[1](inputLayer(X[i]), param)
        X[i+1]  = oneStep(X[i], u)
    end

    return X
end

function oneControllerLoss(x0, param::AbstractArray{T}; totalTimeStep = totalTimeStep) where {T<:Real}
    ltotal  = 0.0f0
    for xi in x0
        X = Vector{Vector{T}}(undef, totalTimeStep+1)
        X[1] = xi
        for i in 1:totalTimeStep
            ui = controlArray[1](inputLayer(X[i]), param)
            X[i+1] = oneStep(X[i], ui)
        end
        ltotal += setDistanceLoss(X) 
    end
    return ltotal/length(x0)
end

function testBayesian(xi, param; totalTimeStep = totalTimeStep)    
    X =  trajectory(xi, param; totalTimeStep = totalTimeStep)
    clf()
    plots(X, fig1)
    # plotPartition(X, ψ, θk)
    return X
end

function trainOneController()
    param       = 0.1f0*randn(Float32, controlNN_length[1])
    opt         = Adam(0.001f0)
    counter     = 0
    minibatch   = 3
    diff_results = DiffResults.GradientResult(param)

    for i in 1:10000

        x0      = sampleInitialState(param; totalTimeStep=8000, minibatch=minibatch)
        l1(θ)   = oneControllerLoss(x0, θ; totalTimeStep = 1500)

        ForwardDiff.gradient!(diff_results, l1, param)
        grads = DiffResults.gradient(diff_results)

        if counter > 5
            if rand() > 0.8
                xi = [0.0f0, pi, 0.0f0, 0.2f0]
            else
                xi = x0[1]
            end
            X = testBayesian(xi, param; totalTimeStep=10000)
            println("loss = ", oneControllerLoss([xi], param))
            counter = 0
        end
        counter += 1
        Flux.update!(opt, param, grads)
    end
end
