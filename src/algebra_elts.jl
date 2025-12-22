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

function promote_object end

maybe_promote(a, _, ::Nothing) = a, nothing
maybe_promote(a, b, map) = promote_with_map(a, b, map)

function promote_with_map(a::StarAlgebra, mstr, map)
    new_obj = promote_object(object(a), mstr, map)
    return StarAlgebra(new_obj, mstr), map
end

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
