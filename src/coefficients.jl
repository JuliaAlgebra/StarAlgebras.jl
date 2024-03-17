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

"""
    canonical(ac::AbstractCoefficients)
Compute the canonical form of `ac` (e.g. grouping coefficients together, etc).

If `ac` can be brough to canonical form in-place one has to implement
* `MA.mutability(::Type{typeof(ac)}, canonical, ::Vararg{Type}) = MA.IsMutable()`
* `MA.operate!(canonical, ac)` that performs this canonicalization.

otherwise `canonical(ac)` needs to be implemented.
"""
canonical(ac::AbstractCoefficients) = ac
MA.operate(::typeof(canonical), x) = canonical(x) # fallback?

# example implementation for vectors
function MA.mutability(
    ::Type{<:SparseVector},
    ::typeof(canonical),
    ::Vararg{Type},
)
    return MA.IsMutable()
end
MA.operate!(::typeof(canonical), sv::SparseVector) = dropzeros!(sv)
canonical(v::AbstractVector) = v

function Base.:(==)(ac1::AbstractCoefficients, ac2::AbstractCoefficients)
    ac1 = MA.operate!!(canonical, ac1)
    ac2 = MA.operate!!(canonical, ac2)
    return keys(ac1) == keys(ac2) && values(ac1) == values(ac2)
end

"""
    nonzero_pairs(ac::AbstractCoefficients)
Return an iterator over pairs `(k=>v)` of keys and values stored in `ac`.

The iterator contains all pairs with `v` potentially non-zero.
"""
function nonzero_pairs(ac::AbstractCoefficients)
    return (k => v for (k, v) in zip(keys(ac), values(ac)))
end

nonzero_pairs(v::AbstractVector) = pairs(v)
function nonzero_pairs(v::AbstractSparseVector)
    return zip(SparseArrays.nonzeroinds(v), SparseArrays.nonzeros(v))
end

aug(ac::AbstractCoefficients) = sum(c * aug(x) for (x, c) in pairs(ac))
aug(v::AbstractVector) = sum(v)
aug(x::Any) = 1 # ???? dubious...
