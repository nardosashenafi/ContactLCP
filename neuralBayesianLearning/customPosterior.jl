using Random
using AdvancedVI
using Flux
struct PosteriorDistribution{T} <: ContinuousMultivariateDistribution
    mψ::Vector{T}
    σψ::Vector{T}
    θk::Vector{T}

    function PosteriorDistribution(mψ::Vector{T}, σψ::Vector{T}, θk::Vector{T}) where {T<:Real}
        mψ = mψ
        σψ = σψ
        θk = θk
        new{T}(mψ, σψ, θk)
    end
end

function Distributions.logpdf(d::PosteriorDistribution, y::Vector{T}) where {T<:Real}
    sum(vcat(logpdf.(Normal.(d.mψ, d.σψ), y[1:length(d.mψ)]), 
        logpdf.(Bernoulli.(Flux.sigmoid(d.θk)), y[length(d.mψ)+1:end])))
end

function Distributions.rand(rng::AbstractRNG, d::PosteriorDistribution)
    rand.(vcat(Normal.(d.mψ, d.σψ), Bernoulli.(Flux.sigmoid(d.θk))))
end

function Distributions.length(d::PosteriorDistribution)
    return sum(length(mψ) + length(θk))
end

function Distributions.entropy(d::PosteriorDistribution, z)
    -logpdf(d, z)
end

function Base.max(d::PosteriorDistribution)

    bernoulliMax(θk) = θk >= 0.5
    return vcat(d.mψ, bernoulliMax.(θk))
end

function (elbo::ELBO)(
    rng::AbstractRNG,
    alg::ADVI,
    q::AdvancedVI.VariationalPosterior,
    logπ::Function,
    num_samples
)
    #   𝔼_q(z)[log p(xᵢ, z)]
    # = ∫ log p(xᵢ, z) q(z) dz
    # = ∫ log p(xᵢ, f(ϕ)) q(f(ϕ)) |det J_f(ϕ)| dϕ   (since change of variables)
    # = ∫ log p(xᵢ, f(ϕ)) q̃(ϕ) dϕ                   (since q(f(ϕ)) |det J_f(ϕ)| = q̃(ϕ))
    # = 𝔼_q̃(ϕ)[log p(xᵢ, z)]

    #   𝔼_q(z)[log q(z)]
    # = ∫ q(f(ϕ)) log (q(f(ϕ))) |det J_f(ϕ)| dϕ     (since q(f(ϕ)) |det J_f(ϕ)| = q̃(ϕ))
    # = 𝔼_q̃(ϕ) [log q(f(ϕ))]
    # = 𝔼_q̃(ϕ) [log q̃(ϕ) - log |det J_f(ϕ)|]
    # = 𝔼_q̃(ϕ) [log q̃(ϕ)] - 𝔼_q̃(ϕ) [log |det J_f(ϕ)|]
    # = - ℍ(q̃(ϕ)) - 𝔼_q̃(ϕ) [log |det J_f(ϕ)|]

    # Finally, the ELBO is given by
    # ELBO = 𝔼_q(z)[log p(xᵢ, z)] - 𝔼_q(z)[log q(z)]
    #      = 𝔼_q̃(ϕ)[log p(xᵢ, z)] + 𝔼_q̃(ϕ) [log |det J_f(ϕ)|] + ℍ(q̃(ϕ))

    # If f: supp(p(z | x)) → ℝ then
    # ELBO = 𝔼[log p(x, z) - log q(z)]
    #      = 𝔼[log p(x, f⁻¹(z̃)) + logabsdet(J(f⁻¹(z̃)))] + ℍ(q̃(z̃))
    #      = 𝔼[log p(x, z) - logabsdetjac(J(f(z)))] + ℍ(q̃(z̃))

    # But our `forward(q)` is using f⁻¹: ℝ → supp(p(z | x)) going forward → `+ logjac`
    _, z, logjac, _ = forward(rng, q)
    res = (logπ(z) + logjac) / num_samples
    res += entropy(q, z) / num_samples

    for i = 2:num_samples
        _, z, logjac, _ = forward(rng, q)
        res += (logπ(z) + logjac) / num_samples
        res += entropy(q, z) / num_samples
    end

    return res
end
