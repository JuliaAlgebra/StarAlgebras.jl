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

"""
    ExplicitBasis
Explicit bases are stored in an `AbstractVector` and hence immutable
(e.g. fixed in length).
"""
abstract type ExplicitBasis{T,I} <: AbstractBasis{T,I} end

function coeffs(
    cfs::AbstractCoefficients,
    source::AbstractBasis,
    target::AbstractBasis,
)
    source === target && return cfs
    res = SparseCoefficients(eltype(target)[], valtype(cfs)[])
    return coeffs!(res, cfs, source, target)
end
