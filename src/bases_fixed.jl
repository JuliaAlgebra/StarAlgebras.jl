# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

abstract type FiniteSupportBasis{T,I} <: ExplicitBasis{T,I} end

"""
    supp(fb::FiniteSupportBasis)
Return the supporting elements of `fb` as an indexable vector
"""
function supp end
supp(fb::FiniteSupportBasis) = fb.supporting_elts

Base.IteratorSize(::Type{<:FiniteSupportBasis}) = Base.HasLength()
Base.IteratorEltype(::Type{<:FiniteSupportBasis{T}}) where {T} = T
Base.length(b::FiniteSupportBasis) = length(supp(b))

Base.iterate(b::FiniteSupportBasis) = iterate(supp(b))
Base.iterate(b::FiniteSupportBasis, state) = iterate(supp(b), state)
# function Base.IndexStyle(::Type{<:FiniteSupportBasis{T,I,V}}) where {T,I,V}
#     return Base.IndexStyle(V)
# end

Base.@propagate_inbounds function Base.getindex(
    b::FiniteSupportBasis,
    i::Integer,
)
    return supp(b)[i]
end

mutable struct FixedBasis{T,I,V<:AbstractVector{T}} <:
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
Base.in(x, b::FixedBasis) = haskey(b.relts, x)
Base.getindex(b::FixedBasis{T}, x::T) where {T} = b.relts[x]
Base.getindex(b::FixedBasis, i::Integer) = b.supporting_elts[i]

struct SubBasis{T,I,K,V<:AbstractVector{K},B<:AbstractBasis{T,K}} <:
       FiniteSupportBasis{T,I}
    supporting_idcs::V
    parent_basis::B
    function SubBasis(supporting_idcs::AbstractVector{K}, parent_basis::AbstractBasis{T,K}) where {T,K}
        return new{T,keytype(supporting_idcs),K,typeof(supporting_idcs),typeof(parent_basis)}(supporting_idcs, parent_basis)
    end
end

supp(sb::SubBasis) = sb.supporting_idcs
Base.parent(sub::SubBasis) = sub.parent_basis

Base.in(x, b::SubBasis) = b.parent_basis[x] in supp(b)
function Base.getindex(b::SubBasis{T,I}, x::T) where {T,I}
    return convert(I, findfirst(isequal(b.parent_basis[x]), supp(b)))
end