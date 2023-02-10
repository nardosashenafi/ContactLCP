using LinearAlgebra, Flux, DiffEqFlux, LogExpFunctions, ForwardDiff
using Distributions
using PyPlot, BSON, LaTeXStrings, FileIO, MATLAB 
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
vis = startAnimator()
include("trainingHelpers.jl")
include("testHelpers.jl")

const tspan         = (0.0f0, 4.0f0) 
const Δt            = 0.001f0
const totalTimeStep = Int(floor(tspan[2]/Δt))
const binSize       = 3

binNN = FastChain(FastDense(5, 4, elu),
                    FastDense(4, 3, elu),
                    FastDense(3, binSize))

const binNN_length = DiffEqFlux.paramlength(binNN) 

controlArray     = Array{Function}(undef, binSize);
controlNN_length = Vector{Int}(undef, binSize)

for i in 1:binSize 
    controlArray[i] = FastChain(FastDense(5, 10, elu),
                            FastDense(10, 4, elu),
                            FastDense(4, 1))
    
    controlNN_length[i] = DiffEqFlux.paramlength(controlArray[i]) 
end

function softargmax(x; β=50.0f0) 
    return exp.(β*x)/sum(exp.(β*x))
end

function bin(x, ψ::AbstractArray{T}) where {T<:Real} 
    return Flux.softmax(binNN(inputLayer(x), ψ))
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

function setDistancelossPerState(x)
    x1, x2, x1dot, x2dot = x

    return 15.0f0*(1.0f0-cos(x2)) + 
            0.8f0*abs(x1dot) + 0.8f0*abs(x2dot)
end

function setDistanceLoss(X::Vector{Vector{T}}, pk, k; r=1e-10) where {T<:Real}
    lmin, lminIndex = findmin(map(setDistancelossPerState, X))

    actionProbability = 1.0f0 - mean(map((pki, ki) -> pki[ki], pk[1:lminIndex], k[1:lminIndex]))
    actionProbability > 0.05 ? delta = lmin*actionProbability : delta = lmin

    incurCostAt = (TRACK_LENGTH)/2.0 
    doubleHingeLoss = 0.0f0
    for x in X 
        abs(x[1]) > incurCostAt ? doubleHingeLoss += 2.0f0*(abs(x[1])-incurCostAt) : nothing
    end

    # return ifelse(delta < r, 0.0f0, delta - r) + doubleHingeLoss 
    return delta + doubleHingeLoss 

end

function averageControlLoss(x0, param::AbstractArray{T}, explorationPercent; totalTimeStep = totalTimeStep) where {T<:Real}
    ψ, θk   = unstackParams(param)
    ltotal = 0.0f0

    for xi in x0
        X = Vector{Vector{T}}(undef, totalTimeStep+1)
        X[1] = xi

        pk = Vector{Vector{T}}(undef, totalTimeStep)
        k = Vector{Int}(undef, totalTimeStep)

        for i in 1:totalTimeStep
            pk[i] = bin(X[i], ψ)      #Gaussian: quadratic boudaries

            if rand() > explorationPercent
                k[i] = argmax(pk[i])
            else
                k[i] = rand(1:1:binSize)
            end

            ui = input(X[i], θk, k[i])
            X[i+1] = oneStep(X[i], ui)
        end
        ltotal += setDistanceLoss(X[1:totalTimeStep], pk, k) 
       end
    return ltotal/length(x0)
end

function poi(x, ψ)
    return bin(x, ψ)
end

function trainEM()

    ψ           = 0.5f0*randn(Float64, binNN_length)
    θk          = [0.1f0*randn(Float64, controlNN_length[i]) for i in 1:binSize]
    param       = stackParams(ψ, θk)
    opt         = Adam(0.001f0)
    counter     = 0
    minibatch   = 3
    diff_results = DiffResults.GradientResult(param)
    initialExploration = 0.60

    for i in 1:90000
        ψ, θk   = unstackParams(param)
        #Apply combination of DAgger and uniform sampling around the desired equilbrium state
        #inorder to assist exploration 
        x0      = sampleInitialState(ψ, θk; totalTimeStep=8000, minibatch=minibatch)

        # For each state in the trajectory, compute the loss incurred by each of the given by the 
        # bin generator 
        explorationPercent = exp.(-i*0.0005)*initialExploration + 0.05
        l1(θ)   = averageControlLoss(x0, θ, explorationPercent; totalTimeStep=1200)
        # l1(θ)   = accumulatedLossSampler(x0, θ, explorationPercent; totalTimeStep=1500)

        # ForwardDiff.gradient takes the gradient of the loss wrt the state bin parameters and the 
        # controller parameters in each bin
        # lg1     = ForwardDiff.gradient(l1, param)

        ForwardDiff.gradient!(diff_results, l1, param)
        grads = DiffResults.gradient(diff_results)

        if counter > 10
            ψ, θk  = unstackParams(param)
            if rand() > 0.8
                xi = [0.0f0, pi, 0.0f0, 0.2f0]
            else
                xi = x0[1]
            end
            X = testBayesian(xi, ψ, θk; totalTimeStep=7000)
            println("loss = ", averageControlLoss([xi], param, explorationPercent; totalTimeStep=7000), " POI = ", poi(xi, ψ))
            # BSON.@save "neuralBL/savedWeights/setdistancetraining.bson" param
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