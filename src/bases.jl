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
expressed in it.
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

function zero_coeffs(::Type{S}, ::ImplicitBasis{T,I}) where {S,T,I}
    return SparseCoefficients(I[], S[])
end

"""
    ExplicitBasis{T,I}
Explicit bases are stored e.g. in an `AbstractVector` and hence immutable
(in particular of well defined and fixed length).
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
function coeffs(cfs, source::AbstractBasis, target::AbstractBasis)
    source === target && return cfs
    source == target && return cfs
    res = zero_coeffs(value_type(cfs), target)
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
