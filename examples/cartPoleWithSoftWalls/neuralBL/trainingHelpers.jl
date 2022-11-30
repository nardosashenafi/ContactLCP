
oneStep(x, param) = ContactLCP.oneTimeStep(lcp, x, param; Δt=Δt)[1]

function integrate(xi, ψ::Vector{T}, θk; totalTimeStep = totalTimeStep) where {T<:Real}
    X = Vector{Vector{T}}()
    U = Vector{T}()

    pk = Vector(undef, totalTimeStep)
    c = Vector{Int}(undef, totalTimeStep)

    push!(X, xi)

    for i in 1:totalTimeStep
        pk[i] = bin(X[end], ψ)           #holds a softmax
        c[i]  = argmax(pk[i])             # can be categorical
        u     = input(X[end], θk, c[i])
        x     = oneStep(X[end], u)
        push!(X, x)
        push!(U, u[1])
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
    
    x0 = T.([0.0, rand(-pi:0.005:pi), 
         rand(-1.0:0.001:1.0), rand(-1.0:0.001:1.0)])


    X  = trajectory(x0, ψ, θk)
    samples = Vector{Vector{T}}()
    for i in 1:minibatch
        if rand() > 0.2
            select      = rand(X)
            select[1]   = rand(-0.5f0:0.01f0:0.5f0)   #stop the training from visiting large values of x
            push!(samples, select)
        else
            push!(samples, [rand(-0.1:0.001:0.1), rand(-0.2:0.001:0.2), 
                            rand(-0.5:0.001:0.5), rand(-0.5:0.001:0.5)])
        end
    end

    return samples
end

