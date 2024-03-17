"""
    abstract type AbstractCoefficients{V,K} end
Implements `Base.keys`, `Base.values`.
"""
abstract type AbstractCoefficients{K,V} end

key_type(::Type{<:AbstractCoefficients{K}}) where {K} = K
Base.valtype(::Type{<:AbstractCoefficients{K,V}}) where {K,V} = V
key_type(b::AbstractCoefficients) = key_type(typeof(b))
Base.valtype(b::AbstractCoefficients) = valtype(typeof(b))

key_type(b) = keytype(b)
# `keytype(::Type{SparseVector{V,K}})` is not defined so it falls
# back to `keytype{::Type{<:AbstractArray})` which returns `Int`.
key_type(::Type{SparseArrays.SparseVector{V,K}}) where {V,K} = K
key_type(v::SparseArrays.SparseVector) = key_type(typeof(v))

Base.iszero(ac::AbstractCoefficients) = isempty(keys(ac))

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
