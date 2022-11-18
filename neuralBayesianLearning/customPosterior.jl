using Random
using AdvancedVI
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
    # yn = y[1:length(d.mψ)]
    # yb = y[length(d.mψ)+1:end]

    # return sum(vcat(logpdf(MvNormal(d.mψ, d.σψ), yn), 
    #         logpdf.(Bernoulli.(Flux.sigmoid.(d.θk)), yb)))
    return sum(logpdf.(Bernoulli.(berAct.(d.θk)), y))
end

function Distributions.rand(rng::AbstractRNG, d::PosteriorDistribution)
    # rn = rand(MvNormal(d.mψ, Flux.softplus.(d.σψ)))
    rb = Float32.(rand.(Bernoulli.(berAct.(d.θk))))

    # return vcat(rn, rb)
end

function Distributions.length(d::PosteriorDistribution)
    # return sum(length(d.mψ) + length(d.θk))
    return length(d.θk)
end

function Distributions.entropy(d::PosteriorDistribution)
    # entropy(MvNormal(d.mψ, Flux.softplus.(d.σψ))) + sum([entropy(Bernoulli(Flux.sigmoid.(θ))) for θ in d.θk])
    return sum([entropy(Bernoulli(berAct.(θ))) for θ in d.θk])
end

function Base.max(d::PosteriorDistribution)

    bernoulliMax(θk) = θk >= 0.5
    # return vcat(d.mψ, bernoulliMax.(Flux.sigmoid.(d.θk)))
    return bernoulliMax.(berAct.(d.θk))
end


