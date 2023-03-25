include("EMSwingup.jl")

function sampleInitialState(param::Vector{T}; totalTimeStep = totalTimeStep, minibatch=4) where {T<:Real}
    
    x0  = T.([rand(-TRACK_LENGTH/2.0:0.01:TRACK_LENGTH/2.0), rand(0.0:0.005:2pi), 
                  rand(-0.5:0.001:0.5), rand(-3.0:0.001:3.0)])

    X           = trajectory(x0, param)
    samples     = Vector{Vector{T}}(undef, minibatch)
    cropBuffer  = 0.3f0

        for i in 1:minibatch
        if rand() > 0.4
            select     = rand(X)
            #It is highly likely that the initial parameter cause the cart to go far from x1 = 0
            #Hence sampling from this trajectory tends to overwhelm the reply buffer with large x1 values.
            #Clamping x1 is bad because it tends to always work under the edge cases
            #so simply randomize x1 when it goes too far
            if !(-TRACK_LENGTH/2.0+w/2.0 <= select[1] <= TRACK_LENGTH/2.0-w/2.0)
                select[1]  = rand(-TRACK_LENGTH/2.0+w/2.0:0.01:TRACK_LENGTH/2.0-w/2.0) 
            end
            samples[i] = select
        else
            samples[i] =  [rand(d:0.01:D+d), rand(-0.5:0.001:0.5), 
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

function setDistanceLoss(X::Vector{Vector{T}}; r=0.1f0) where {T<:Real}
    delta, _ = findmin(map(setDistancelossPerState, X))

    incurCostAt = (TRACK_LENGTH)/2.0 
    doubleHingeLoss = 0.0f0
    for x in X 
        abs(x[1]) > incurCostAt ? doubleHingeLoss += 2.0f0*(abs(x[1])-incurCostAt) : nothing
    end

    return delta + doubleHingeLoss
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
    # animate(X)
    # sleep(1)
    plotPartition(X, param)
    return X
end

function trainOneController()
    param       = 0.3f0*randn(Float32, controlNN_length[1])
    opt         = Adam(0.001f0)
    counter     = 0
    minibatch   = 3
    diff_results = DiffResults.GradientResult(param)

    for i in 1:10000

        x0      = sampleInitialState(param; totalTimeStep=8000, minibatch=minibatch)
        l1(θ)   = oneControllerLoss(x0, θ; totalTimeStep = 1200)

        ForwardDiff.gradient!(diff_results, l1, param)
        grads = DiffResults.gradient(diff_results)

        if counter > 20
            if rand() > 0.8
                xi = [0.0f0, pi, 0.0f0, 0.2f0]
            else
                xi = x0[1]
            end
            X = testBayesian(xi, param; totalTimeStep=10000)
            println("One controller loss = ", oneControllerLoss([xi], param))
            counter = 0
        end
        counter += 1
        Flux.update!(opt, param, grads)
    end
end


function plotPartition(X, param)
    width = 30
    
    X2    = range(-2.0f0pi, 2.0f0pi, length=width)
    X2dot = range(-6.0f0, 6.0f0, length=width)

    u = Matrix{Float32}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(D+d, X2[j], 0.0f0, X2dot[width-i+1])
            u[i, j] = controlArray[1](inputLayer(x), param)[1]
        end
    end

    plt.subplot(2, 2, 2)
    imshow(u, extent = [X2[1], X2[end], X2dot[1], X2dot[end]])
    ylabel("Control")
    xlabel("x2 vs x2dot")

    #overlap the trajectory on the control heat map
    plot(getindex.(X, 2), getindex.(X, 4), "r")
    scatter(X[end][2], X[end][4])


end