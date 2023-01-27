include("EMSwingup.jl")

function oneBatch(xi, param::AbstractArray{T}; totalTimeStep = totalTimeStep) where {T<:Real}
    
    ψ, θk   = unstackParams(param)
    ltotal  = 0.0f0

    X = Vector{Vector{T}}(undef, totalTimeStep+1)
    X[1] = xi
    for i in 1:totalTimeStep
        pk = bin(X[i], ψ)    

        uall = [input(X[i], θk, k)[1] for k in 1:binSize]
        ui = [dot(pk, uall)]
        X[i+1] = oneStep(X[i], ui)
    end
    ltotal += setDistanceLoss(X) + 0.1f0*lossPerState(X[end]) 
    return ltotal
end

function computeLoss(x0, param::AbstractArray{T}; totalTimeStep = totalTimeStep) where {T<:Real}
    
    l = Threads.Atomic{T}()
    Threads.@threads for xi in x0
        Threads.atomic_add!(l, oneBatch(xi, param))
    end
    return l.value/length(x0)
end

function gradient!(grad, x0, param::AbstractArray{T}; totalTimeStep = 1000) where {T<:Real}
    Threads.@threads for i in eachindex(x0)
        grad[i] = ForwardDiff.gradient((θ) -> oneBatch(x0[i], θ; totalTimeStep = totalTimeStep), param)
    end
end

function trainParallel()

    ψ           = 0.5f0*randn(Float32, binNN_length)
    θk          = [0.1f0*randn(Float32, controlNN_length[i]) for i in 1:binSize]
    param       = stackParams(ψ, θk)
    opt         = Adam(0.001f0)
    counter     = 0
    minibatch   = 4

    for i in 1:10000
        ψ, θk   = unstackParams(param)
        x0      = sampleInitialState(ψ, θk; totalTimeStep=8000, minibatch=minibatch)

        lg1     = Vector{Vector{eltype(param)}}(undef, length(x0))
        gradient!(lg1, x0, param; totalTimeStep = 1500) 
        grads   = mean(lg1)

        if counter > 10
            ψ, θk  = unstackParams(param)
            if rand() > 0.8
                xi = [0.0f0, pi, 0.0f0, 0.2f0]
            else
                xi = deepcopy(x0[1])
            end
            testBayesian(xi, ψ, θk; totalTimeStep=7000)
            l1   = oneBatch(xi, param; totalTimeStep = 7000)
            println("loss = ", l1, " POI = ", poi(xi, ψ))
            BSON.@save "neuralBL/savedWeights/beastHangingWalls.bson" param
            counter = 0
        end
        counter += 1
        Flux.update!(opt, param, grads)
    end
end