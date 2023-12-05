"""
    AbstractBasis{T,I}
Implements a bijection between basis elements and integers.

 * `in(x, basis)` for basis elements
 * `Base.getindex(A, i::I) → T`
 * `Base.getindex(A, t::T) → I` # the bijection part

"""
abstract type AbstractBasis{T,I} end

Base.eltype(::Type{<:AbstractBasis{T}}) where {T} = T
Base.keytype(::Type{<:AbstractBasis{T,I}}) where {T,I} = I
Base.keytype(b::AbstractBasis) = SparseArrays.indtype(typeof(b))

"""
    ImplicitBasis{T,I}
Implicit bases are not stored in memory and can be potentially infinite.
"""
abstract type ImplicitBasis{T,I} <: AbstractBasis{T,I} end

"""
    ExplicitBasis
Explicit bases are stored in an AbstractVector and hence immutable
(e.g. fixed in length)
"""
abstract type ExplicitBasis{T,I} <: AbstractBasis{T,I} end

struct Basis{T,I,A<:AbstractVector{T}} <: ExplicitBasis{T,I}
    basis::A
    rbasis::Dict{T,I}
end

function Basis{I}(basis::AbstractVector) where {I}
    Base.require_one_based_indexing(basis)
    length(basis) <= typemax(I) ||
        throw("index type $I is to small for basis of length $(length(basis))")
    @assert !(eltype(basis) <: Integer)
    return Basis(basis, Dict(b => I(idx) for (idx, b) in pairs(basis)))
end

Base.size(b::Basis) = size(b.basis)
Base.IndexStyle(::Type{<:Basis{T,I,A}}) where {T,I,A} = Base.IndexStyle(A)

Base.@propagate_inbounds Base.getindex(b::Basis{T,I}, i::I) = b.basis[i]
Base.@propagate_inbounds Base.getindex(b::Basis{T}, g::T) where {T} = b.rbasis[g]

Base.in(g, b::Basis) = haskey(b.rbasis, g)

# convenience only:
Base.@propagate_inbounds function Base.getindex(b::Basis, i::Integer)
    idx = convert(SparseArrays.indtype(b), i)
    return b[idx]
end
