# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

"""
    mutable struct FixedBasis{T,I,V<:AbstractVector{T}} <: ExplicitBasis{T,I}
        supporting_elts::V
        relts::Dict{T,I}
        starof::Vector{I}
    end
"""
mutable struct FixedBasis{T,I,V<:AbstractVector{T}} <: ExplicitBasis{T,I}
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
Base.@propagate_inbounds Base.getindex(b::FixedBasis, i::Integer) = b.supporting_elts[i]

Base.iterate(b::FixedBasis) = iterate(b.supporting_elts)
Base.iterate(b::FixedBasis, state) = iterate(b.supporting_elts, state)
function Base.IndexStyle(::Type{<:FixedBasis{T,I,V}}) where {T,I,V}
    return Base.IndexStyle(V)
end

struct SubBasis{T,I,K,V<:AbstractVector{K},B<:AbstractBasis{T,K}} <:
       ExplicitBasis{T,I}
    supporting_idcs::V
    parent_basis::B
    is_sorted::Bool
    function SubBasis(supporting_idcs::AbstractVector{K}, parent_basis::AbstractBasis{T,K}) where {T,K}
        return new{T,keytype(supporting_idcs),K,typeof(supporting_idcs),typeof(parent_basis)}(supporting_idcs, parent_basis, issorted(supporting_idcs))
    end
end

Base.parent(sub::SubBasis) = sub.parent_basis

Base.length(b::SubBasis) = length(b.supporting_idcs)
function _iterate(b::SubBasis, elem_state)
    if isnothing(elem_state)
        return
    end
    elem, state = elem_state
    return parent(b)[elem], state
end
Base.iterate(b::SubBasis) = _iterate(b, iterate(b.supporting_idcs))
Base.iterate(b::SubBasis, st) = _iterate(b, iterate(b.supporting_idcs, st))

function Base.get(b::SubBasis{T,I}, x::T, default) where {T,I}
    key = b.parent_basis[x]
    if b.is_sorted
        i = searchsortedfirst(b.supporting_idcs, key)
        if i in eachindex(b.supporting_idcs) && b.supporting_idcs[i] == key
            return convert(I, i)
        end
    else
        i = findfirst(isequal(key), b.supporting_idcs)
        if !isnothing(i)
            return convert(I, i)
        end
    end
    return default
end

Base.in(x::T, b::SubBasis{T}) where T = !isnothing(get(b, x, nothing))

Base.getindex(b::SubBasis, i::Integer) = parent(b)[b.supporting_idcs[i]]
function Base.getindex(b::SubBasis{T,I}, x::T) where {T,I}
    i = get(b, x, nothing)
    if isnothing(i)
        throw(BoundsError(b, x))
    end
    return i::I
end