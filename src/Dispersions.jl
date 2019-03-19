module Dispersions

export AbstractDispersion,
TwoLevelSystem

abstract type AbstractDispersion end

struct TwoLevelSystem <: AbstractDispersion
    D₀::Float64
    k₀::Float64
    γp::Float64
    ω₀::Float64

    TwoLevelSystem(D₀=0., k₀=Inf, γp=1e8, ω₀=k₀) = new(float(D₀), float(k₀), float(γp), float(ω₀))
end

(tls::TwoLevelSystem)(k::Number) = tls.γp/(k-tls.k₀+1im*tls.γp)
function (tls::TwoLevelSystem)(k::Number,ks::Array,ψ::Array)
    γk = tls(k)
    h = zeros(Float64,size(ψ,1))
    for i ∈ eachindex(ks)
        h += abs2(ks[i])*abs2.(ψ[:,i])
    end
    return γk./(1 .+ h)
end

end # module
