# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

Base.hash(a::AlgebraElement, h::UInt) = hash(coeffs(a), hash(parent(a), h))

function Base.:(==)(X::AlgebraElement, Y::AlgebraElement)
    parent(X) == parent(Y) || return false
    return coeffs(X) == coeffs(Y)
end

Base.copy(a::AlgebraElement) = AlgebraElement(copy(coeffs(a)), parent(a))

function Base.deepcopy_internal(a::AlgebraElement, id::IdDict)
    if !haskey(id, a)
        id[a] = AlgebraElement(Base.deepcopy_internal(coeffs(a), id), parent(a))
    end
    return id[a]
end

# call overload:
(a::AlgebraElement)(x) = coeffs(a)[basis(a)[x]]
Base.setindex!(a::AlgebraElement, v, idx) = a.coeffs[basis(a)[idx]] = v

# AlgebraElement specific functions

function supp(a::AlgebraElement)
    b = basis(a)
    return [b[i] for (i, _) in nonzero_pairs(coeffs(a))]
end

function LinearAlgebra.norm(a::AlgebraElement, p::Real=2)
    return LinearAlgebra.norm(coeffs(a), p)
end

function LinearAlgebra.dot(a::AlgebraElement, b::AlgebraElement)
    return LinearAlgebra.dot(coeffs(a), coeffs(b))
end

function LinearAlgebra.dot(w::AbstractVector, b::AlgebraElement)
    return LinearAlgebra.dot(w, coeffs(b))
end
function LinearAlgebra.dot(a::AlgebraElement, w::AbstractVector)
    @assert key_type(basis(parent(a))) <: Integer
    return LinearAlgebra.dot(coeffs(a), w)
end

"""
    promote_object(object, mstructure, map)

Promote the `object` field of a `StarAlgebra` to the new
multiplicative structure `mstructure` and the key map `map`.
"""
function promote_object end

"""
    maybe_promote(a, parent, map)

Promote `a` given the promoted parent `parent` and the key map
`map. Call this function as follows:

    maybe_promote(a, ::Nothing)

Just return `a` in case there is no promotion to do.

To be used as follows
```julia
function promote_basis_with_maps(a::..., b::...)
    _a, _b = promote_basis_with_maps(parent(a), parent(b))
    return maybe_promote(a, _a...), maybe_promote(b, _b...)
end
```
"""
function maybe_promote end

maybe_promote(a, _, ::Nothing) = a, nothing
maybe_promote(a, b, map) = promote_with_map(a, b, map)

"""
    promote_with_map(a, parent, map)

Promote `a` given the promoted parent `parent` and the key map
`map. Implement this function but call `maybe_promote`.
"""
function promote_with_map end

function promote_with_map(a::StarAlgebra, mstr, map)
    new_obj = promote_object(object(a), mstr, map)
    return StarAlgebra(new_obj, mstr), map
end

"""
    promote_basis_with_map(a, b)

Return `(new_a, map_a), (new_b, map_b)` where `new_a` and `new_b`
are the promoted version of `a` and `b` so that they now have the same basis.
The function `map_a` (resp. `map_b`) maps the key of `a` (resp. `b`) to the
corresponding key of `new_a` (resp. `new_b`).

The convention for implementing this interface for a type `A` with parent type `B`
and function `parent` that obtains the parent type is as follows:
```julia
import StarAlgebras as SA

function SA.promote_with_map(a::A, b::B, map)
    # Promote other fields of `a` using `b` and `map`...
    return A(b, ...)
end

function SA.promote_basis_with_maps(a::A, b::A)
    _a, _b = SA.promote_basis_with_maps(parent(a), parent(b))
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end

function SA.promote_basis_with_maps(a::A, b::B)
    _a, _b = SA.promote_basis_with_maps(parent(a), b)
    return SA.maybe_promote(a, _a...), _b
end
```
"""
function promote_basis_with_maps end

function promote_basis_with_maps(
    a::StarAlgebra,
    b::StarAlgebra,
)
    _a, _b = promote_basis_with_maps(mstructure(a), mstructure(b))
    return maybe_promote(a, _a...), maybe_promote(b, _b...)
end

function promote_basis_with_maps(
    a::StarAlgebra,
    b::Union{MultiplicativeStructure,AbstractBasis},
)
    _a, _b = promote_basis_with_maps(mstructure(a), b)
    return maybe_promote(a, _a...), _b
end

function promote_with_map(a::AlgebraElement, alg, map)
    c = coeffs(a)
    if map !== identity
        c = map_keys(map, c)
    end
    return AlgebraElement(c, alg), map
end

function promote_basis_with_maps(
    a::AlgebraElement,
    b::AlgebraElement,
)
    _a, _b = promote_basis_with_maps(parent(a), parent(b))
    return maybe_promote(a, _a...), maybe_promote(b, _b...)
end

function promote_basis_with_maps(
    a::AlgebraElement,
    b::Union{
        StarAlgebra,
        MultiplicativeStructure,
        AbstractBasis,
    },
)
    _a, _b = promote_basis_with_maps(parent(a), b)
    return maybe_promote(a, _a...), _b
end

"""
    promote_basis(a, b)

Same as `promote_basis_with_maps` but without the key maps.
"""
function promote_basis(a, b)
    _a, _b = promote_basis_with_maps(a, b)
    return _a[1], _b[1]
end
