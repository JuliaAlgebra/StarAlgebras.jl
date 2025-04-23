# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

"""
    mutable struct DiracBasis{T,F} <: ImplicitBasis{T,T}
        object::S
    end

Given any iterable `object` over elements of type `T`, the basis
corresponds to the list of elements of `object`.
The type `S` should implement the methods:
```julia
Base.IteratorSize(::Type{S})
Base.eltype(::S) # Should return `T`
Base.iterate(::S)
Base.iterate(::S, state)
Base.in(::T, ::S)
```
Optionally, it may implement the following,
```julia
Base.isempty(::S)
Base.length(::S)
```
"""
struct DiracBasis{T,S} <: ImplicitBasis{T,T}
    object::S

    function DiracBasis(itr)
        @assert !isempty(itr)
        return new{eltype(itr),typeof(itr)}(itr)
    end
end

function Base.IteratorSize(::Type{<:DiracBasis{T,S}}) where {T,S}
    return Base.IteratorSize(S)
end

Base.iterate(db::DiracBasis) = iterate(object(db))
Base.iterate(db::DiracBasis, st) = iterate(object(db), st)

function Base.getindex(db::DiracBasis{T}, x::T) where {T}
    @assert x in object(db)
    return x
end

"""
    mutable struct MappedBasis{T,I,S,F,U} <: ImplicitBasis{T,I}
        object::S
        map::F
        inverse_map::U
    end

Given any iterable `object` over elements of type `I`, the basis
corresponds to the list of elements `map(i)` corresponding
to the indices `i` of `object`.
The type `S` should implement the methods:
```julia
Base.IteratorSize(::Type{S})
Base.eltype(::S) # Should return `I`
Base.iterate(::S)
Base.iterate(::S, state)
Base.in(::I, ::S)
```
Optionally, it may implement the following,
```julia
Base.isempty(::S)
Base.length(::S)
```
"""
mutable struct MappedBasis{T,I,S,F,U} <: ImplicitBasis{T,I}
    object::S
    map::F
    inverse_map::U

    function MappedBasis{T,I}(itr, map, inverse_map) where {T,I}
        @assert !isempty(itr)
        @assert T != I || (map == inverse_map == identity)
        return new{T,I,typeof(itr),typeof(map),typeof(inverse_map)}(itr, map, inverse_map)
    end
end

function Base.IteratorSize(::Type{<:MappedBasis{T,I,S}}) where {T,I,S}
    return Base.IteratorSize(S)
end

function MappedBasis(itr, map, inverse_map)
    return MappedBasis{typeof(map(first(itr))),eltype(itr)}(itr, map, inverse_map)
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
    return b.inverse_map(x)
end

function Base.getindex(b::MappedBasis{T,I}, x::I) where {T,I}
    return b.map(x)
end

object(db::Union{DiracBasis,MappedBasis}) = db.object

function Base.size(b::Union{DiracBasis,MappedBasis})
    @assert Base.haslength(object(b))
    return size(object(b))
end

function Base.length(b::Union{DiracBasis,MappedBasis})
    @assert Base.haslength(object(b))
    return length(object(b))
end
