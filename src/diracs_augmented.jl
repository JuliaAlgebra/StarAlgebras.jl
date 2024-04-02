struct Augmented{K} <: AbstractCoefficients{K,Int}
    elt::K
end # corresponds to (elt - 1)

function Base.getindex(aδ::Augmented{K}, i::K) where {K}
    w = aug(aδ.elt)
    isone(aδ.elt) && return zero(w)
    i == aδ.elt && return w
    isone(i) && return -w
    return zero(w)
end

canonical(aδ::Augmented) = aδ
Base.iszero(aδ::Augmented) = all(iszero, values(aδ))

Base.keys(aδ::Augmented) = (k = keys(aδ.elt); (one(first(k)), first(k)))
function Base.values(aδ::Augmented)
    (e, x) = keys(aδ)
    return (aδ[e], aδ[x])
end

function Base.show(io::IO, aδ::Augmented)
    ioc = IOContext(io, :limit => true)
    if iszero(aδ)
        print(ioc, "(0)")
    else
        e, x = keys(aδ)
        _, v = values(aδ)
        print(ioc, '(', -v, '·', e, '+', v, '·', x, ')')
    end
end

Base.isless(ad1::Augmented, ad2::Augmented) = isless(ad1.elt, ad2.elt)

aug(::Augmented) = 0
