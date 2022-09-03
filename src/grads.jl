
include("lcp.jl")

function ChainRulesCore.frule((_, _), 
    ::typeof(lcpSolve), lcp::Lcp, x::Vector{T}, θ::Vector{T}) where {T<:Real}

    model = JuMP.Model(() -> DiffOpt.diff_optimizer(Mosek.Optimizer))
    sol, A, b = lcpSolve(lcp, x, θ, model=model)

    Δcons1, Δcons2, Δcons3, Δobj = getPerturbations(lcp, model, x, θ, A, b, sol)
    cons1 = model[:cons1]
    cons2 = model[:cons2]
    cons3 = model[:cons3]

    if !isequal(Δcons1, nothing)
        MOI.set.(  
            model, 
            DiffOpt.ForwardConstraintFunction(),   
            cons1, 
            Δcons1
        )
    end

    if !isequal(Δcons2, nothing)
        MOI.set.(  
            model, 
            DiffOpt.ForwardConstraintFunction(),    
            cons2, 
            Δcons2
        )
    end

    if !isequal(Δcons3, nothing)
        MOI.set.(  
            model, 
            DiffOpt.ForwardConstraintFunction(),  
            cons3, 
            Δcons3
        )
    end

    if !isequal(Δobj, nothing)
        MOI.set.(  
            model, 
            DiffOpt.ForwardObjectiveFunction(),    
            Δobj
        )
    end

    DiffOpt.forward_differentiate!(JuMP.backend(model))

    ∂sol∂θm =  MOI.get.(   
                model,
                DiffOpt.ForwardVariablePrimal(),
                [model[:q]; model[:v]; model[:λ]])

    return (sol, ∂sol∂θm)
end

function getPerturbations(lcp::Lcp, model, x, θ, A, b, sol)
    #provide the derivative of the objective and constraints of the LCP wrt the learned parameters.
    #this tells DiffOpt which parameter we are interested in differentiating the solution with respect to.
    return lcp.sys(model, x, θ, A, b, sol)   
end

function lossGrad(lcp::Lcp, S, Λ, θm::Vector{T}, loss::Function, ∂l∂optsol::Function, Qs::Matrix{T}, Qλ::Matrix{T}) where {T<:Real}

    l   = 0.0
    lg  = 0.0

    for i in 1:lcp.sys.totalTimeStep
        sol1, grad1     = ChainRulesCore.frule((nothing, nothing), lcpSolve, lcp, S[i], θm)
        s, λ            = unstackSol(lcp, sol1)
        l               += loss(S[i], s, Qs, Λ[i], λ, Qλ)
        lg              += ∂l∂optsol(S[i], s, Qs, Λ[i], λ, Qλ) * grad1
    end

    # loss(θm) = computeLoss(lcp, S, Λ, θm, x0, stateNum; Δt = Δt, totalTimeStep=totalTimeStep)
    # grad_m = ForwardDiff.derivative(loss, θm0)

    return l, lg
end

function checkGradient(lcp::Lcp, x1::Vector{T}, θ1::Vector{T}, θ2::Vector{T}) where {T<:Real}

    sol1, grad1 =  ChainRulesCore.frule((nothing, nothing), lcpSolve, lcp, x1, θ1)
    q1, v1, λ1 = (sol1[1:stateNum], sol1[stateNum+1:2stateNum], sol1[2stateNum+1:end])

    sol2, grad2 =  ChainRulesCore.frule((nothing, nothing), lcpSolve, lcp, x1, θ2)
    q2, v2, λ2 = (sol2[1:stateNum], sol2[stateNum+1:2stateNum], sol2[2stateNum+1:end])

    finiteDiff_q = (q2 - q1) / (θ2 - θ1)
    error_q = abs.(finiteDiff_q - grad1[1:stateNum])

    finiteDiff_v = (v2 - v1) / (θ2 - θ1)
    error_v = abs.(finiteDiff_v - grad1[stateNum+1:2stateNum])

    finiteDiff_λ = (λ2 - λ1) / (θ2 - θ1)
    error_λ = abs.(finiteDiff_λ - grad1[2stateNum+1:end])

    return error_q, error_v, error_λ, finiteDiff_q, finiteDiff_v, finiteDiff_λ, grad1
    # return error_q./finiteDiff_q, error_v./finiteDiff_v, grad1
end