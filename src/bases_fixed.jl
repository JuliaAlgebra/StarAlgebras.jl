# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

mutable struct FixedBasis{T,I,V<:AbstractVector{T}} <:
               ExplicitBasis{T,I}
    elts::V
    relts::Dict{T,I}
    starof::Vector{I}
end

"""
    FixedBasis(elts::AbstractVector)

Represents a linear basis which is fixed of length (given by `elts`).
The elements of `elts` are assumed to be _linearly independent_ and (as a set) invariant under `star`.
Due to its finiteness, `b::FixedBasis` can be indexed like a vector (by positive integers) and therefore elements expressed in `b` can be represented by simple (sparse) Vectors.

## Examples

```julia
julia> import StarAlgebras as SA;

julia> SA.star(s::String) = reverse(s) # identity is also fine

julia> b = SA.FixedBasis(["", "a", "b", "aa", "ab", "ba", "bb"]);

julia> b["a"]
1

julia> b["ab"]
4

julia> b[4]
"ab"

julia> b["abc"]
ERROR: KeyError: key "abc" not found
[...]

julia> A = SA.StarAlgebra(String, SA.DiracMStructure(b, *))
*-algebra of String

julia> c = A("a") + A("b")
1·"a" + 1·"b"

julia> c*A("a")
1·"aa" + 1·"ba"

julia> A("a")*c
1·"aa" + 1·"ab"

julia> c*c
1·"aa" + 1·"ab" + 1·"ba" + 1·"bb"

julia> SA.coeffs(c*c)
7-element SparseArrays.SparseVector{Int64, Int64} with 4 stored entries:
  [4]  =  1
  [5]  =  1
  [6]  =  1
  [7]  =  1
```
"""
FixedBasis(elts::AbstractVector{T}) where {T} = FixedBasis{T,keytype(elts)}(elts)

function FixedBasis{T,I}(elts::AbstractVector{T}) where {T,I}
    relts = Dict(b => I(idx) for (idx, b) in pairs(elts))
    starof = [relts[star(x)] for x in elts]
    return FixedBasis{T,I,typeof(elts)}(elts, relts, starof)
end

function FixedBasis{T,I}(basis::AbstractBasis{T}; n::Integer) where {T,I}
    return FixedBasis{T,I}(collect(Iterators.take(basis, n)))
end

FixedBasis(basis::AbstractBasis{T}; n::Integer) where {T} = FixedBasis{T,typeof(n)}(basis; n)

Base.in(x, b::FixedBasis) = haskey(b.relts, x)
Base.getindex(b::FixedBasis{T}, x::T) where {T} = b.relts[x]
Base.getindex(b::FixedBasis, i::Integer) = b.elts[i]

Base.length(b::FixedBasis) = length(b.elts)
Base.iterate(b::FixedBasis) = iterate(b.elts)
Base.iterate(b::FixedBasis, state) = iterate(b.elts, state)
Base.IndexStyle(::Type{<:FixedBasis{T,I,V}}) where {T,I,V} = Base.IndexStyle(V)

"""
    struct SubBasis{T,I,K,B<:AbstractBasis{T,K},V<:AbstractVector{K}} <:
        ExplicitBasis{T,I}
        parent_basis::B
        keys::V
        is_sorted::Bool
    end

Represents a sub-basis of a given basis, where `keys` is a vector of keys
representing the sub-basis elements, `parent_basis` is the parent basis from
which the sub-basis is derived, and `is_sorted` indicates whether the keys are
sorted.
"""
struct SubBasis{T,I,K,B<:AbstractBasis{T,K},V<:AbstractVector{K}} <:
       ExplicitBasis{T,I}
    parent_basis::B
    keys::V
    is_sorted::Bool
    function SubBasis(parent_basis::AbstractBasis{T,K}, keys::AbstractVector{K}) where {T,K}
        return new{T,keytype(keys),K,typeof(parent_basis),typeof(keys)}(parent_basis, keys, issorted(keys))
    end
end

Base.copy(b::SubBasis) = SubBasis(copy(parent(b)), copy(b.keys))

Base.parent(sub::SubBasis) = sub.parent_basis

Base.length(b::SubBasis) = length(b.keys)
function _iterate(b::SubBasis, elem_state)
    if isnothing(elem_state)
        return
    end
    elem, state = elem_state
    return parent(b)[elem], state
end
Base.iterate(b::SubBasis) = _iterate(b, iterate(b.keys))
Base.iterate(b::SubBasis, st) = _iterate(b, iterate(b.keys, st))

function Base.get(b::SubBasis{T,I}, x::T, default) where {T,I}
    key = b.parent_basis[x]
    if b.is_sorted
        i = searchsortedfirst(b.keys, key)
        if i in eachindex(b.keys) && b.keys[i] == key
            return convert(I, i)
        end
    else
        i = findfirst(isequal(key), b.keys)
        if !isnothing(i)
            return convert(I, i)
        end
    end
    return default
end

Base.in(x::T, b::SubBasis{T}) where T = !isnothing(get(b, x, nothing))
Base.haskey(b::SubBasis, i::Integer) = i in eachindex(b.keys)

Base.getindex(b::SubBasis, i::Integer) = parent(b)[b.keys[i]]
function Base.getindex(b::SubBasis{T,I}, x::T) where {T,I}
    i = get(b, x, nothing)
    if isnothing(i)
        throw(KeyError(x))
    end
    return i::I
end
