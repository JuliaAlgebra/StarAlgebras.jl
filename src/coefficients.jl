"""
    abstract type AbstractCoefficients{V,K} end
Implements `Base.keys`, `Base.values`.
"""
abstract type AbstractCoefficients{K,V} end

Base.keytype(::Type{<:AbstractCoefficients{K}}) where {K} = K
Base.valtype(::Type{<:AbstractCoefficients{K,V}}) where {K,V} = V
Base.keytype(b::AbstractCoefficients) = keytype(typeof(b))
Base.valtype(b::AbstractCoefficients) = valtype(typeof(b))

Base.iszero(sc::AbstractCoefficients) = isempty(keys(sc))

function Base.:(==)(ac1::AbstractCoefficients, ac2::AbstractCoefficients)
    for x in (ac1, ac2)
        if !__iscanonical(x)
            __canonicalize!(x)
        end
    end
    return keys(ac1) == keys(ac2) && values(ac1) == values(ac2)
end

aug(ac::AbstractCoefficients) = sum(c * aug(x) for (x, c) in pairs(ac))
aug(v::AbstractVector) = sum(v)
aug(x::Any) = 1 # ???? dubious...
