# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

"""
    mutable struct MappedBasis{T,I,S,F,U} <: ImplicitBasis{T,I}
        object::S
        map::F
        inverse_map::U
    end

Given any iterable `object` over elements of type `I`, the basis
corresponds to the list of elements `map(i)` corresponding
to the indices `i` of `object`.
"""
mutable struct MappedBasis{T,I,S,F,U} <: ImplicitBasis{T,I}
    object::S
    map::F
    inverse_map::U

    function MappedBasis{T,I}(itr, map, inverse_map) where {T,I}
        @assert !isempty(itr)
        return new{T,I,typeof(itr),typeof(map),typeof(inverse_map)}(itr, map, inverse_map)
    end
end

function MappedBasis(itr, map, inverse_map)
    return MappedBasis{typeof(map(first(itr))),eltype(itr)}(itr, map, inverse_map)
end

identity_basis(itr) = MappedBasis(itr, identity, identity)

object(db::MappedBasis) = db.object

function Base.IteratorSize(::Type{<:MappedBasis{T,I,S}}) where {T,I,S}
    return Base.IteratorSize(S)
end
function Base.length(b::MappedBasis)
    @assert Base.haslength(object(b))
    return length(object(b))
end
function _iterate(b::MappedBasis, elem_state)
    if isnothing(elem_state)
        return
    end
    elem, state = elem_state
    return b.map(elem), state
end
Base.iterate(b::MappedBasis) = _iterate(b, iterate(object(b)))
Base.iterate(b::MappedBasis, st) = _iterate(b, iterate(object(b), st))

Base.in(g::T, b::MappedBasis{T}) where {T} = haskey(b, b.inverse_map(g))
Base.haskey(b::MappedBasis{T,I}, k::I) where {T,I} = k in object(b)

function Base.getindex(b::MappedBasis{T}, x::T) where {T}
    return b.unmap(x)
end

function Base.getindex(b::MappedBasis{T,I}, x::I) where {T,I}
    return b.map(x)
end

function Base.getindex(::MappedBasis{T,T,S,typeof(identity),typeof(identity)}, x::T) where {T,S}
    return x
end
