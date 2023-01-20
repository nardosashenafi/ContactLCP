

function integrate(xi, ψ::Vector{T}, θk; totalTimeStep = totalTimeStep) where {T<:Real}
    X = Vector{Vector{T}}(undef, totalTimeStep)
    U = Vector{Vector{T}}(undef, totalTimeStep-1)

    X[1] = xi

    for i in 1:totalTimeStep-1
        pk      = bin(X[i], ψ)           #holds a softmax
        c       = argmax(pk)             # can be categorical
        U[i]    = input(X[i], θk, c)
        X[i+1]  = oneStep(X[i], U[i])
    end

    return X, U
end

function trajectory(x0, ψ, θk; totalTimeStep = totalTimeStep)
    X,_    = integrate(x0, ψ, θk; totalTimeStep = totalTimeStep)
    return X
end

function accumulatedLoss(traj)
    l = 1.0f0/length(traj)*sum(map(x -> dot(x, x), traj))
    return l
end

function sampleInitialState(ψ::Vector{T}, θk; totalTimeStep = totalTimeStep, minibatch=4) where {T<:Real}
    
    X = Vector{Vector{T}}()

    for i in 1:2     
        x0  = T.([rand(d+w/2:0.01:D+d-w/2), rand(0.0:0.005:2pi), 
                rand(-0.5:0.001:0.5), rand(-3.0:0.001:3.0)])
        append!(X, trajectory(x0, ψ, θk))
    end

    samples     = Vector{Vector{T}}(undef, minibatch)
    cropBuffer  = 0.3f0

    for i in 1:minibatch
        if rand() > 0.2
            select     = rand(X)
            #It is highly likely that the initial parameter cause the cart to go far from x1 = 0
            #Hence sampling from this trajectory tends to overwhelm the reply buffer with large x1 values.
            #Clamping x1 is bad because it tends to always work under the edge cases
            #so simply randomize x1 when it goes too far
            if !(d+w/2 - cropBuffer <= select[1] <= D+d-w/2 + cropBuffer)
                select[1]  = rand(d+w/2 - cropBuffer:0.01:D+d-w/2 + cropBuffer) 
            end
            samples[i] = select
        else
            samples[i] =  [rand(-0.1:0.001:0.1), rand(-0.5:0.001:0.5), 
                            rand(-0.5:0.001:0.5), rand(-1.0:0.001:1.0)]
        end
    end

    return samples
end

