# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

"""
    mutable struct DiracBasis{T,F} <: ImplicitBasis{T,T}
        object::S
    end

Given any iterable `object` over elements of type `T`, the basis
corresponds to the list of elements of `object`.
The `object` (its type `S`, respectively) should implement the methods:
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

object(db::DiracBasis) = db.object

function Base.IteratorSize(::Type{<:DiracBasis{T,S}}) where {T,S}
    return Base.IteratorSize(S)
end

function Base.size(db::DiracBasis)
    # Some object such as `PermGroup` do not define `size`
    return (length(db),)
end

function Base.length(db::DiracBasis)
    @assert Base.haslength(object(db))
    return length(object(db))
end

Base.iterate(db::DiracBasis) = iterate(object(db))
Base.iterate(db::DiracBasis, st) = iterate(object(db), st)

Base.in(g, db::DiracBasis) = g in object(db)
Base.haskey(db::DiracBasis, g) = in(g, db)

Base.firstindex(db::DiracBasis) = Base.first(object(db))
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
corresponds to the list of elements `map(i::I)` where `i` 
belongs to `object`. The function `inverse_map` should
be the inverse of `map`. That is, `inverse_map(map(i))` should be `i`
for all `i in object`.
The `object` (its type `S`, respectively) should implement the methods:
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

    function MappedBasis{T}(itr, map, inverse_map) where {T}
        @assert !isempty(itr)
        I = eltype(itr)
        @assert T != I || (map == inverse_map == identity)
        return new{T,I,typeof(itr),typeof(map),typeof(inverse_map)}(itr, map, inverse_map)
    end
end

function MappedBasis(itr, map, inverse_map)
    return MappedBasis{typeof(map(first(itr)))}(itr, map, inverse_map)
end

function Base.:(==)(a::MappedBasis, b::MappedBasis)
    return object(a) == object(b) && a.map == b.map && a.inverse_map == b.inverse_map
end

object(db::MappedBasis) = db.object

function Base.IteratorSize(::Type{<:MappedBasis{T,I,S}}) where {T,I,S}
    return Base.IteratorSize(S)
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

function Base.size(b::MappedBasis)
    @assert Base.haslength(object(b))
    return size(object(b))
end

function Base.length(b::MappedBasis)
    @assert Base.haslength(object(b))
    return length(object(b))
end

Base.in(g::T, b::MappedBasis{T}) where {T} = haskey(b, b.inverse_map(g))
Base.haskey(b::MappedBasis{T,I}, k::I) where {T,I} = k in object(b)

function Base.getindex(b::MappedBasis{T}, x::T) where {T}
    return b.inverse_map(x)
end

function Base.getindex(b::MappedBasis{T,I}, x::I) where {T,I}
    return b.map(x)
end

function promote_basis_with_maps(a::ImplicitBasis, b::ImplicitBasis)
    if a == b
        return (a, nothing), (b, nothing)
    end
    error("Bases $a and $b are different and do not support promotion.")
end
