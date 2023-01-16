using LinearAlgebra, Flux, DiffEqFlux, LogExpFunctions, ForwardDiff
using Distributions
using PyPlot 
using ForwardDiff: Chunk, GradientConfig
# import ContactLCP

using JuMP
include("../../../src/lcp.jl")
include("../../../src/solver.jl")
# include("../hangingWallDynamics.jl")
include("../hangingWallDynamics.jl")

function checkContact(gn::Vector{T}, gThreshold, total_contact_num) where {T<:Real}
     
    whichInContact = zeros(T, total_contact_num)

    for i in 1:total_contact_num
        if (-wallThickness/2.0f0+gThreshold < gn[i] < gThreshold)  
            # we need to detect contact on both sides of a wall. If we simply use gn < gThreshold,
            # then the pendulum would stick to the wall due to the opposing contact forces from both sides
            # of the wall. To combat that, we need to deactivate the contact force from one side of the wall
            # after a certain depth.
            whichInContact[i] = 1 
        end
    end
    current_contact_num = Int(sum(whichInContact))

    return findall(x -> x == 1, whichInContact), current_contact_num
end

# include("../longWallDynamics.jl")

sys  = CartPoleWithSoftWalls()
lcp  = Lcp(Float32, sys)
oneStep(x, param) = oneTimeStep(lcp, x, param; Δt=Δt)

fig1 = figure()
include("trainingHelpers.jl")
include("testHelpers.jl")

const tspan         = (0.0f0, 4.0f0) 
const Δt            = 0.001f0
const totalTimeStep = Int(floor(tspan[2]/Δt))
const binSize       = 3

binNN = FastChain(FastDense(5, 5, elu),
                    FastDense(5, 4, elu),
                    FastDense(4, binSize))

const binNN_length = DiffEqFlux.paramlength(binNN) 

controlArray     = Array{Function}(undef, binSize);
controlNN_length = Vector{Int}(undef, binSize)

for i in 1:binSize 
    controlArray[i] = FastChain(FastDense(5, 8, elu),
                            FastDense(8, 4, elu),
                            FastDense(4, 1))
    
    controlNN_length[i] = DiffEqFlux.paramlength(controlArray[i]) 
end

println("Check the following properties in the initialization")
println("Make sure that bin(x, ψ) does not return a large probability in some state.
        This causes argmax(bin(x, ψ)) to return the same value for a long time in the training
        until the highest probability crosses 0.5")

function bin(x, θ::Vector{T}) where {T<:Real} 
    return Flux.softmax(binNN(inputLayer(x), θ))
end

function input(x, θk, i::Int)
    return controlArray[i](inputLayer(x), θk[i])
end

function unstackControlParam(param::Vector{T}) where {T<:Real}

    θk = Vector{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}}(undef, binSize)
    θk[1] = @view param[1:controlNN_length[1]] 

    for i in 2:binSize
        θk[i] = @view param[sum([controlNN_length[j] for j in 1:i-1])+1:sum([controlNN_length[j] for j in 1:i])] 
    end

    return θk
end

function unstackParams(param)
    return param[1:binNN_length], unstackControlParam(param[binNN_length+1:end])
end

function stackParams(ψ::Vector{T}, θk::Vector{Vector{T}}) where {T<:Real}
    param = Vector{T}()
    param = vcat(param, ψ)
    for i in 1:binSize
        append!(param, θk[i])
    end
    return param
end

function lossPerState(x)
    x1, x2, x1dot, x2dot = x
    doubleHinge_x = 0.0f0

    abs(x1) > D/2.0f0 ? doubleHinge_x = 2.0f0*abs.(x1) : nothing

    # high cost on x1dot to lower fast impact
    return doubleHinge_x  + 12.0f0*(1.0f0-cos(x2)) + 
            2.0f0*x1dot^2.0f0 + 0.1f0*x2dot^2.0f0

end

function computeLossSampler(x0, param::Vector{T}; totalTimeStep = totalTimeStep) where {T<:Real}
    ψ, θk   = unstackParams(param)
    ltotal  = 0.0f0

    for xi in x0
        x = xi
        for i in 1:totalTimeStep
            pk = bin(x, ψ)

            if rand() > 0.3
                k = argmax(pk)      # improve the greedy
            else
                k = rand(1:1:binSize)    # exploration
            end

            x  = oneStep(x, input(x, θk, k) )
            lk = pk[k]*lossPerState(x) 
            ltotal += lk/totalTimeStep
        end
    end
    return ltotal/length(x0)
end

function computeLoss(X, param::Vector{T}, sampleEvery::Int, fullTraj::Bool ; totalTimeStep = totalTimeStep) where {T<:Real}
    ψ, θk   = unstackParams(param)
    ltotal  = 0.0f0

    for i in eachindex(X)
        X2      = X[i][1:sampleEvery:end]
        len     = length(X2)

        for x in X2
            pk     = bin(x, ψ)
            lk     = 0.0f0
            for k in 1:binSize 
                u  = input(x, θk, k)          
                x2 = oneStep(x, u)
                lk += pk[k]*lossPerState(x2) 
            end
            ltotal += lk/len
        end
    end
    return ltotal/length(X)
end

function oneBatch(xi, param::Vector{T}; totalTimeStep = totalTimeStep) where {T<:Real}
    
    ψ, θk   = unstackParams(param)
    ltotal  = 0.0f0
    pk      = Vector(undef, totalTimeStep)
    x2      = deepcopy(xi)

    for i in 1:totalTimeStep
        pk[i]  = bin(x2, ψ)           #holds a softmax
        lk     = 0.0f0
        for k in 1:binSize     
            u  = input(x2, θk, k)          
            x  = oneStep(x2, u)
            lk += pk[i][k]*lossPerState(x)       #maximize the changes of selecting lower loss
        end
        ltotal += lk/totalTimeStep

        ### select the next state
        c  = rand(Categorical(bin(x2, ψ)))
        # c  = argmax(bin(x2, ψ))
        x2 = oneStep(x2, input(x2, θk, c))
    end
    return ltotal
end

function computeLoss(x0, param::Vector{T}; totalTimeStep = totalTimeStep) where {T<:Real}
    
    l = Threads.Atomic{T}()
    Threads.@threads for xi in x0
        Threads.atomic_add!(l, oneBatch(xi, param))
    end
    return 3.0f0*l.value/length(x0)
end

function gradient!(grad, x0, param::Vector{T}; totalTimeStep = 1000) where {T<:Real}
    Threads.@threads for i in eachindex(x0)
        grad[i] = ForwardDiff.gradient((θ) -> oneBatch(x0[i], θ; totalTimeStep = totalTimeStep), param)
    end
end

function poi(ψ)
    xi = [0.0f0, -4.0f0, 0.0f0, 2.0f0]
    return bin(xi, ψ)
end

function trainEM()

    ψ           = 0.5f0*randn(Float32, binNN_length)
    θk          = [0.1f0*randn(Float32, controlNN_length[i]) for i in 1:binSize]
    param       = stackParams(ψ, θk)
    opt         = Adam(0.001f0)
    counter     = 0
    minibatch   = 3
    diff_results = DiffResults.GradientResult(param)

    for i in 1:10000
        ψ, θk   = unstackParams(param)
        #Apply combination of DAgger and uniform sampling around the desired equilbrium state
        #inorder to assist exploration 
        x0      = sampleInitialState(ψ, θk; totalTimeStep=8000, minibatch=minibatch)

        # X       = [trajectory(xi, ψ, θk; totalTimeStep = 3000) for xi in x0]
        # For each state in the trajectory, compute the loss incurred by each of the given by the 
        # bin generator 
        # l1(θ)   = computeLoss(X, θ, 2, true; totalTimeStep = 2000)
        l1(θ)   = computeLossSampler(x0, θ; totalTimeStep = 4000)

        # ForwardDiff.gradient takes the gradient of the loss wrt the state bin parameters and the 
        # controller parameters in each bin
        # lg1     = ForwardDiff.gradient(l1, param)

        ForwardDiff.gradient!(diff_results, l1, param)
        ∇ = DiffResults.gradient(diff_results)
        Δ = Flux.Optimise.apply!(opt, param, ∇)
        @. param = param - Δ

        if counter > 5
            ψ, θk  = unstackParams(param)
            if rand() > 0.8
                xi = [0.0f0, pi, 1.0f0, 0.5f0]
            else
                xi = x0[1]
            end
            X = testBayesian(xi, ψ, θk; totalTimeStep=10000)
            println("loss = ", computeLossSampler([xi], param), " POI = ", poi(ψ))
            counter = 0
        end
        counter += 1
        Flux.update!(opt, param, lg1)
    end
end

function trainParallel()

    ψ           = 0.5f0*randn(Float32, binNN_length)
    θk          = deepcopy(ps)
    param       = stackParams(ψ, θk)
    opt         = Adam(0.001f0)
    counter     = 0
    minibatch   = 6

    for i in 1:10000
        ψ, θk   = unstackParams(param)
        x0      = sampleInitialState(ψ, θk; totalTimeStep=5000, minibatch=minibatch)

        lg1     = Vector{Vector{eltype(param)}}(undef, length(x0))
        gradient!(lg1, x0, param; totalTimeStep = 1000) 
        grads   = mean(lg1)

        if counter > 5
            ψ, θk  = unstackParams(param)
            if rand() > 0.8
                xi = [0.0f0, pi, 1.0f0, 0.5f0]
            else
                xi = deepcopy(x0[1])
            end
            testBayesian(xi, ψ, θk; totalTimeStep=7000)
            l1   = oneBatch(xi, param; totalTimeStep = 7000)
            println("loss = ", l1)
            counter = 0
        end
        counter += 1
        Flux.update!(opt, param, grads)
    end
end

function checkGradient(param)
    x0  = [[0.0f0, pi, 1.0f0, 0.5f0]]

    l1(θ)  = computeLoss(x0, θ; totalTimeStep = 1000)
    lg1    = ForwardDiff.gradient(l1, param)

    δ = 0.001
    fd = zeros(length(param))

    for ind in 1:length(param)
        param2 = deepcopy(param)
        param2[ind] += δ

        fd[ind] = (l1(param2) - l1(param))/δ
    end

    percError = abs.(lg1 - fd) ./ fd

    println("error = ", percError)
end