using Random
using AdvancedVI, Distributions
using Flux

struct PosteriorDistribution{T} <: ContinuousMultivariateDistribution
    θk::Vector{T}

    function PosteriorDistribution(θk::Vector{T}) where {T<:Real}
        θk = θk
        new{T}(θk)
    end
end

function berAct(θk)
    return Flux.sigmoid(θk)
end

function Distributions.logpdf(d::PosteriorDistribution, y::Vector{T}) where {T<:Real}

    return sum(logpdf.(Bernoulli.(berAct.(d.θk)), y))
end

function Distributions.rand(rng::AbstractRNG, d::PosteriorDistribution)
    rb = Float32.(rand.(Bernoulli.(berAct.(d.θk))))
end

function Distributions.length(d::PosteriorDistribution)
    return length(d.θk)
end

function Distributions.entropy(d::PosteriorDistribution)
    return sum([entropy(Bernoulli(berAct.(θ))) for θ in d.θk])
end

function Bijectors.bijector(d::PosteriorDistribution)
    return Logit.(0.0f0, 1.0f0)     #for bernoully ∂z/∂θk ≈ 1
end

function Base.max(d::PosteriorDistribution)

    if d.θk >= 0.5f0
        return 0
    else 
        return 1
    end

end


