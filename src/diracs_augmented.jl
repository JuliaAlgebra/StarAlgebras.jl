struct AugmentedDirac{K,V} <: AbstractCoefficients{K,V}
    dirac::Dirac{K,V}
end
AugmentedDirac(x) = AugmentedDirac(Dirac(x, Int(!isone(x))))

function Base.getindex(aδ::AugmentedDirac{K}, i::K) where {K}
    v = aδ.dirac[i]
    v = ifelse(isone(aδ.dirac.element), zero(v), v)
    return ifelse(isone(i), -v, v)
end

__iscanonical(aδ::AugmentedDirac) = true
Base.iszero(aδ::AugmentedDirac) = all(iszero, values(aδ))

Base.keys(aδ::AugmentedDirac) = (k = keys(aδ.dirac); (one(first(k)), first(k)))
function Base.values(aδ::AugmentedDirac)
    (e, x) = keys(aδ)
    return (aδ[e], aδ[x])
end

function Base.show(io::IO, aδ::AugmentedDirac)
    ioc = IOContext(io, :limit => true)
    e, x = keys(aδ)
    _, v = values(aδ)
    print(ioc, '(', -v, '·', e, '+', v, '·', x, ')')
end

Base.isless(ad1::AugmentedDirac, ad2::AugmentedDirac) = isless(ad1.dirac, ad2.dirac)

aug(::AugmentedDirac) = 0
