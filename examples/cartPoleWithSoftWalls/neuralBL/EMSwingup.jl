using LinearAlgebra, Flux, DiffEqFlux, LogExpFunctions, ForwardDiff
using Distributions
using PyPlot, BSON
using ForwardDiff: Chunk, GradientConfig
# import ContactLCP

using JuMP
include("../../../src/lcp.jl")
include("../../../src/solver.jl")
include("../hangingWallDynamics.jl")

function checkContact(x, gn::Vector{T}, gThreshold, total_contact_num) where {T<:Real}
     
    whichInContact = zeros(T, 8)

    for i in 1:8
        if (-wallThickness/4.0f0+gThreshold < gn[i] < gThreshold)  
            # we need to detect contact on both sides of a wall. If we simply use gn < gThreshold,
            # then the pendulum would stick to the wall due to the opposing contact forces from both sides
            # of the wall. To combat that, we need to deactivate the contact force from one side of the wall
            # after a certain depth.
            whichInContact[i] = 1 
        end
    end
    contactIndex = findall(x -> x == 1, whichInContact)

    pendulum_xy = pendulumPos(x[1], x[2])

    threshold = 2*gThreshold

    if (lbl[1] + threshold <= pendulum_xy[1] <= lbr[1] - threshold) &&
        (lbl[2] - threshold <= pendulum_xy[2] <= lbl[2] + threshold)
        [filter!(e -> e ≠ j, contactIndex) for j in 1:4]
        append!(contactIndex, 9)
    end

    if (rbl[1] + threshold <= pendulum_xy[1] <= rbr[1] - threshold) &&
        (rbl[2] - threshold <= pendulum_xy[2] <= rbl[2] + threshold)
        [filter!(e -> e ≠ j, contactIndex) for j in 5:8]
        append!(contactIndex, 10)
    end

    if any(i -> i ∈ contactIndex, [3,4])    #the pendulum cannot contact on the edge and along the length at the same time. Or else it will get stick to the wall
        #[3, 4] corresponds to the edge of the link contacting the wall
        #If this contact is active while the contacts [1,2] are activated, then the edge of the link is right underneathh the wall surface
        [filter!(e -> e ≠ j, contactIndex) for j in 1:2]
    end

    if any(i -> i ∈ contactIndex, [7,8])
        [filter!(e -> e ≠ j, contactIndex) for j in 5:6]
    end

    return contactIndex, length(contactIndex)
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

function bin(x, θ::AbstractArray{T}) where {T<:Real} 
    return Flux.softmax(binNN(inputLayer(x), θ))
end

function input(x, θk, i::Int)
    return controlArray[i](inputLayer(x), θk[i])
end

function unstackControlParam(param::AbstractArray{T}) where {T<:Real}

    θk = Vector{SubArray{T, 1, typeof(param), Tuple{UnitRange{Int64}}, true}}(undef, binSize)
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

    abs(x1) > D/2.0f0 ? doubleHinge_x = 3.0f0*abs.(x1) : nothing

    # high cost on x1dot to lower fast impact and to encourage pumping
    return doubleHinge_x  + 12.0f0*(1.0f0-cos(x2)) + 
            2.0f0*x1dot^2.0f0 + 0.4f0*x2dot^2.0f0

end

function computeLossSampler(x0, param::AbstractArray{T}; totalTimeStep = totalTimeStep) where {T<:Real}
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
            # if i == totalTimeStep   #terminal loss
            #     lk *= 2.0f0
            # end
            ltotal += lk/totalTimeStep
        end
    end
    return ltotal/length(x0)
end

function computeLoss(X, param::AbstractArray{T}, sampleEvery::Int, fullTraj::Bool ; totalTimeStep = totalTimeStep) where {T<:Real}
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

function poi(ψ)
    xi = [0.0f0, 0.0f0, 0.0f0, 0.0f0]
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
        grads = DiffResults.gradient(diff_results)

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