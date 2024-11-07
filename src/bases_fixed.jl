# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

abstract type FiniteSupportBasis{T,I} <: ExplicitBasis{T,I} end

"""
    supp(fb::FiniteSupportBasis)
Return the supporting elements of `fb` as an indexable vector
"""
function supp end

Base.IteratorSize(::Type{<:FiniteSupportBasis}) = Base.HasLength()
Base.IteratorEltype(::Type{<:FiniteSupportBasis{T}}) where {T} = T
Base.length(b::FiniteSupportBasis) = length(supp(b))

Base.iterate(b::FiniteSupportBasis) = iterate(supp(b))
Base.iterate(b::FiniteSupportBasis, state) = iterate(supp(b), state)
# function Base.IndexStyle(::Type{<:FiniteSupportBasis{T,I,V}}) where {T,I,V}
#     return Base.IndexStyle(V)
# end

# To break ambiguity
# ??
Base.@propagate_inbounds function Base.getindex(
    b::FiniteSupportBasis,
    i::Integer,
)
    return supp(b)[i]
end

mutable struct FixedBasis{T,I,V<:AbstractVector{T},M<:MTable{T,I}} <:
               FiniteSupportBasis{T,I}
    supporting_elts::V
    relts::Dict{T,I}
    starof::Vector{I}
end

function FixedBasis{T,I}(elts::AbstractVector{T}) where {T,I}
    relts = Dict(b => I(idx) for (idx, b) in pairs(elts))
    starof = [relts[star(x)] for x in elts]
    return FixedBasis{T,I,typeof(elts)}(elts, relts, starof)
end

FixedBasis(elts::AbstractVector{T}) where {T} = FixedBasis{T,keytype(elts)}(elts)

function FixedBasis{T,I}(basis::AbstractBasis{T}; n::Integer) where {T,I}
    return FixedBasis{T,I}(collect(Iterators.take(basis, n)))
end

FixedBasis(basis::AbstractBasis{T}; n::Integer) where {T} = FixedBasis{T,typeof(n)}(basis; n)
supp(fb::FixedBasis) = fb.supporting_elts
Base.in(x, b::FixedBasis) = haskey(b.relts, x)
Base.getindex(b::FixedBasis{T}, x::T) where {T} = b.relts[x]
Base.getindex(b::FixedBasis, i::Integer) = b.elts[i]

struct SubBasis{T,I,V,B<:AbstractBasis{T,I}} <: FiniteSupportBasis{T,I}
    supporting_elts::V
    parent_basis::B
end

supp(sb::SubBasis) = sb.supporting_elts
Base.parent(sub::SubBasis) = sub.parent_basis

Base.in(x, b::SubBasis) = x in supp(b)
function Base.getindex(b::SubBasis{T,I}, x::T) where {T,I<:Integer}
    k = findfirst(==(x), supp(b))
    isnothing(k) && throw("x=$x is not supported on SubBasis")
    @info T I
    return convert(I, k)
end

function Base.getindex(b::SubBasis{T,T}, x::T) where {T}
    x in supp(b) && return x
    throw("x=$x is not supported on SubBasis")
end

struct SubMStructure{SB<:SubBasis,MS} <: MultiplicativeStructure
    basis::SB
    mstructure::MS
end

function mstructure(b::SubBasis)
    pb = b.parent_basis
    return SubMStructure(b, mstructure(pb))
end

function (mstr::SubMStructure)(x::T, y::T) where {T}
    b = mstr.basis
    xy = mstr.mstructure(b[x], b[y])
    return xy
end

# this is used for key-lookup in mstructures.jl
# MA.operate!(op::UnsafeAddMul, …)
_key(mstr::SubMStructure, k) = findfirst(==(k), supp(mstr.basis))
