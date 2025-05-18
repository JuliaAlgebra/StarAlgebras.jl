# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

"""
    AbstractBasis{T,I}
Implements

 * `in(x, basis)` for basis elements
 * `Base.getindex(b::AbstractBasis{T,I}, idx::I)::T`
 * iteration protocol

`AbstractBasis` from the outside may look like a vector with arbitrary indices,
however it may be possibly infinite (and therefore does not define its length).
The elements of the basis are iterated in increasing order according to `key_isless(b)`.

When `I<:Integer` one may store coefficients w.r.t. the basis in an ordinary
vector, however it's length has no relation to the size of the basis (which may
still be infinite).

* [`ImplicitBasis`](@ref) is best suited for bases which can be enumerated but
  are in princle infinite in size (e.g. monomial basis for the polynomial ring).
  An example implementation of [`DiracBasis`](@ref) wrapping an iterator can be
  used for this purpose.
* [`ExplicitBasis`](@ref) can be used for fixed, finite size bases, in particular
  [`FixedBasis`](@ref) implements `AbstractVector`-type storage for elements.
  This can be used e.g. for representing polynomials as (sparse) vectors of
  coefficients w.r.t. a given `FixedBasis`.
"""
abstract type AbstractBasis{T,I} end

Base.eltype(::Type{<:AbstractBasis{T}}) where {T} = T
Base.eltype(b::AbstractBasis) = eltype(typeof(b))
key_type(::Type{<:AbstractBasis{T,I}}) where {T,I} = I
key_type(b::AbstractBasis) = key_type(typeof(b))

"""
    abstract type ImplicitBasis{T,I} <: AbstractBasis{T,I} end

Implicit bases are bases that contains the product of all its elements.
This makes these bases particularly useful to work with [`AlgebraElement`](@ref)s with supports that can not be reasonably bounded.
Note that these bases may not explictly store its elements in memory as they may be potentially infinite.
"""
abstract type ImplicitBasis{T,I} <: AbstractBasis{T,I} end

function zero_coeffs(::Type{S}, ::ImplicitBasis{T,I}) where {S,T,I}
    return SparseCoefficients(I[], S[])
end

"""
    abstract type ExplicitBasis{T,I<:Integer} <: AbstractBasis{T,I} end

Explicit bases are bases of finite length for which the keys are integers.
"""
abstract type ExplicitBasis{T,I<:Integer} <: AbstractBasis{T,I} end

Base.IteratorSize(::Type{<:ExplicitBasis}) = Base.HasLength()
Base.keys(eb::ExplicitBasis) = Base.OneTo(length(eb))
Base.size(eb::ExplicitBasis) = (length(eb),)

function zero_coeffs(::Type{S}, eb::ExplicitBasis{T,I}) where {S,T,I}
    return spzeros(S, I, length(eb))
end

function Base.getindex(eb::ExplicitBasis{T}, range::AbstractRange{<:Integer}) where {T}
    return T[eb[i] for i in range]
end

"""
    coeffs(cfs, source, target)
Translate coefficients `cfs` in `source::AbstractBasis` to basis
`target::AbstractBasis`.
"""
function coeffs(cfs, source::AbstractBasis, target::AbstractBasis)
    source === target && return cfs
    source == target && return cfs
    res = zero_coeffs(value_type(cfs), target)
    return coeffs!(res, cfs, source, target)
end

function coeffs!(res, cfs, source::AbstractBasis, target::AbstractBasis)
    MA.operate!(zero, res)
    lt = key_isless(target)
    for (k, v) in nonzero_pairs(cfs)
        x = source[k]
        kt = target[x]
        vt = getindex_sorted(res, kt; lt) + v
        setindex_sorted!(res, vt, target[x]; lt)
    end
    return res
end

"""
    adjoint_coeffs(cfs, source, target)
Return `A' * cfs` where `A` is the linear map applied by
`coeffs`.
"""
function adjoint_coeffs(cfs, source::AbstractBasis, target::AbstractBasis)
    source === target && return cfs
    source == target && return cfs
    res = zero_coeffs(value_type(cfs), source)
    return adjoint_coeffs!(res, cfs, source, target)
end

function adjoint_coeffs!(res, cfs, source::AbstractBasis, target::AbstractBasis)
    MA.operate!(zero, res)
    for (k, v) in nonzero_pairs(cfs)
        x = target[k]
        # If `x` is not in `source` then the corresponding row in `A` is zero
        # so the column in `A'` is zero hence we can ignore it.
        if x in source
            res[source[x]] += v
        end
    end
    return res
end
