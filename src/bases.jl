"""
    AbstractBasis{T,I}
Implements a bijection between basis elements and integers.

 * `in(x, basis)` for basis elements
 * `Base.getindex(A, i::I) → T`
 * `Base.getindex(A, t::T) → I` # the bijection part

"""
abstract type AbstractBasis{T,I} end

Base.eltype(::Type{<:AbstractBasis{T}}) where {T} = T
Base.eltype(b::AbstractBasis) = eltype(typeof(b))
key_type(::Type{<:AbstractBasis{T,I}}) where {T,I} = I
key_type(b::AbstractBasis) = key_type(typeof(b))

"""
    ImplicitBasis{T,I}
Implicit bases are not stored in memory and can be potentially infinite.
"""
abstract type ImplicitBasis{T,I} <: AbstractBasis{T,I} end

function zero_coeffs(::Type{S}, ::ImplicitBasis{T}) where {S,T}
    return SparseCoefficients(T[], S[])
end

"""
    ExplicitBasis
Explicit bases are stored in an `AbstractVector` and hence immutable
(e.g. fixed in length).
"""
abstract type ExplicitBasis{T,I} <: AbstractBasis{T,I} end

Base.keys(ib::ExplicitBasis) = Base.OneTo(length(ib))

function zero_coeffs(::Type{S}, eb::ExplicitBasis{T,I}) where {S,T,I}
    return spzeros(S, I, length(eb))
end

function Base.getindex(eb::ExplicitBasis, range::AbstractRange{<:Integer})
    return [eb[i] for i in range]
end

"""
    coeffs(cfs, source, target)
Translate coefficients `cfs` in `source::AbstractBasis` to basis
`target::AbstractBasis`.
"""
function coeffs(
    cfs,
    source::AbstractBasis,
    target::AbstractBasis,
)
    source === target && return cfs
    source == target && return cfs
    res = zero_coeffs(valtype(cfs), target)
    return coeffs!(res, cfs, source, target)
end

function coeffs!(res, cfs, source::AbstractBasis, target::AbstractBasis)
    MA.operate!(zero, res)
    for (k, v) in nonzero_pairs(cfs)
        x = source[k]
        res[target[x]] += v
    end
    return res
end
