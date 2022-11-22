
function oneStep(x, ρ::T) where {T<:Real}
    A1 = [0.0f0 -1.0f0; 2.0f0 0.0f0]
    A2 = [0.0f0 -2.0f0; 1.0f0 0.0f0]

    dx = (1.0f0 - ρ)*A1*x + ρ*A2*x
    return x + dx*Δt
end

function trajectory(x0, ψ, uk; totalTimeStep = totalTimeStep)
    X,_    = integrate(x0, ψ, uk; totalTimeStep = totalTimeStep)
    return X
end

function integrate(xi, ψ::Vector{T}, uk; totalTimeStep = totalTimeStep) where {T<:Real}
    X = Vector{Vector{T}}()
    U = Vector{T}()

    pk = Vector(undef, totalTimeStep)
    c = Vector{Int}(undef, totalTimeStep)

    push!(X, xi)

    for i in 1:totalTimeStep
        pk[i] = bin(X[end], ψ)           #holds a softmax
        c[i]  = argmax(pk[i])             # can be categorical
        x     = oneStep(X[end], uk[c[i]])
        push!(X, x)
        push!(U, uk[c[i]])
    end

    return X, U
end

function accumulatedLoss(traj)
    l = 1.0f0/length(traj)*sum(map(x -> dot(x, x), traj))
    return l
end

function sampleInitialState(ψ::Vector{T}, uk; totalTimeStep = totalTimeStep, minibatch=4) where {T<:Real}
    
    x0 = [rand(-4.0:0.001:4.0), rand(-4.0:0.001:4.0)]
    X  = trajectory(x0, ψ, uk)
    
    samples = Vector{Vector{T}}()
    for i in 1:minibatch
        if rand() > 0.1
            push!(samples, rand(X))
        else
            push!(samples, [rand(-0.1:0.001:0.1), rand(-0.1:0.001:0.1)])
        end

    end

    return samples
end

function desiredTrajectory(x0)

    function expert(x)
        if prod(x) <= 0
            p = 0
        else
            p = 1
        end

        return p
    end

    X = Vector{Vector{Float32}}()
    push!(X, x0)

    for i in 1:totalTimeStep-1
        p    = expert(X[end])
        x    = oneStep(X[end], p)
        push!(X, x)
    end
    plt.subplot(2, 2, 1)
    plot(getindex.(X, 1), getindex.(X, 2))
    l = loss(X)

    width = 30
    X1 = range(-5.0f0, 5.0f0, length=width)
    X2 = range(-5.0f0, 5.0f0, length=width)
    
    u = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            p = expert([X1[j], X2[width-i+1]])
            u[i, j] = p
        end
    end

    imshow(u, extent = [X1[1], X1[end], X2[1], X2[end]])
    ylabel("expert")

    return l
end
