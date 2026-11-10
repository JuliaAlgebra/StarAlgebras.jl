# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# An element of an algebra whose basis is a `SubBasis` is equal to the
# corresponding element of the algebra whose basis is its parent so the hash
# is computed from the latter.
function Base.hash(a::AlgebraElement, h::UInt)
    b = implicit(a)
    return hash(coeffs(b), hash(parent(b), h))
end

function Base.:(==)(X::AlgebraElement, Y::AlgebraElement)
    if parent(X) != parent(Y)
        X = implicit(X)
        Y = implicit(Y)
        parent(X) == parent(Y) || return false
    end
    return coeffs(X) == coeffs(Y)
end

Base.copy(a::AlgebraElement) = AlgebraElement(copy(coeffs(a)), parent(a))

function MA.mutability(::Type{<:AlgebraElement{T,A,C}}) where {T,A,C}
    return _coefficients_mutability(C)
end

function MA.mutable_copy(a::AlgebraElement)
    return AlgebraElement(MA.mutable_copy(coeffs(a)), parent(a))
end

function MA.promote_operation(
    ::typeof(zero),
    ::Type{A},
) where {A<:AlgebraElement}
    return similar_type(A, eltype(A))
end

function remove_leading_term(a::AlgebraElement)
    return MA.operate!(remove_leading_term, copy(a))
end

function MA.operate!(::typeof(remove_leading_term), a::AlgebraElement)
    MA.operate!(remove_leading_term, coeffs(a))
    return a
end

function MA.operate(::typeof(remove_leading_term), a::AlgebraElement)
    return remove_leading_term(a)
end

function MA.promote_operation(
    ::typeof(remove_leading_term),
    ::Type{A},
) where {A<:AlgebraElement}
    return A
end

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

function LinearAlgebra.norm(a::AlgebraElement, p::Real = 2)
    return LinearAlgebra.norm(coeffs(a), p)
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
`map`. Call this function as follows:

    maybe_promote(a, ::Nothing)

Just return `a` in case there is no promotion to do.

To be used as follows
```julia
function promote_bases_with_maps(a::..., b::...)
    _a, _b = promote_bases_with_maps(parent(a), parent(b))
    return maybe_promote(a, _a...), maybe_promote(b, _b...)
end
```
"""
function maybe_promote end

maybe_promote(a, _, ::Nothing) = a, nothing
maybe_promote(a, b, map) = promote_with_map(a, b, map)

"""
    promote_parents(a, b, _a, _b)

Same as `(maybe_promote(a, _a...), maybe_promote(b, _b...))` except that `a`
(resp. `b`) is reused as the promoted version of `b` (resp. `a`) if it does not
need to be promoted. This is only valid if `a` and `b` are parents, e.g., bases,
multiplicative structures or algebras, as their promoted versions are then
equal. This is not the case for algebra elements and terms for instance.
"""
function promote_parents(a, b, _a, _b)
    if isnothing(last(_a))
        return (a, nothing), (a, last(_b))
    elseif isnothing(last(_b))
        return (b, last(_a)), (b, nothing)
    end
    return maybe_promote(a, _a...), maybe_promote(b, _b...)
end

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
    promote_bases_with_map(a, b)

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

function SA.promote_bases_with_maps(a::A, b::A)
    _a, _b = SA.promote_bases_with_maps(parent(a), parent(b))
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end

function SA.promote_bases_with_maps(a::A, b::B)
    _a, _b = SA.promote_bases_with_maps(parent(a), b)
    return SA.maybe_promote(a, _a...), _b
end
```
If `A` is itself a parent, e.g., a basis, a multiplicative structure or an
algebra, use [`promote_parents`](@ref) in the method for `(a::A, b::A)`.
"""
function promote_bases_with_maps end

function promote_bases_with_maps(a::StarAlgebra, b::StarAlgebra)
    a == b && return ((a, nothing), (a, nothing))
    _a, _b = promote_bases_with_maps(mstructure(a), mstructure(b))
    return promote_parents(a, b, _a, _b)
end

function promote_bases_with_maps(
    a::StarAlgebra,
    b::Union{MultiplicativeStructure,AbstractBasis},
)
    _a, _b = promote_bases_with_maps(mstructure(a), b)
    return maybe_promote(a, _a...), _b
end

function promote_with_map(a::AlgebraElement, alg, map)
    c = coeffs(a)
    if map !== identity
        c = map_keys(map, c)
    end
    return AlgebraElement(c, alg), map
end

function promote_with_map(t::Term, alg, map)
    return Term(alg, map(t.index), coefficient(t)), map
end

function promote_bases_with_maps(
    a::Union{AlgebraElement,Term},
    b::Union{AlgebraElement,Term},
)
    _a, _b = promote_bases_with_maps(parent(a), parent(b))
    return maybe_promote(a, _a...), maybe_promote(b, _b...)
end

function promote_bases_with_maps(
    a::Union{AlgebraElement,Term},
    b::Union{StarAlgebra,MultiplicativeStructure,AbstractBasis},
)
    _a, _b = promote_bases_with_maps(parent(a), b)
    return maybe_promote(a, _a...), _b
end

"""
    promote_bases(a, b)

Same as `promote_bases_with_maps` but without the key maps.
"""
function promote_bases(a, b)
    _a, _b = promote_bases_with_maps(a, b)
    return _a[1], _b[1]
end

"""
    implicit_basis(basis::AbstractBasis)

Return the implicit basis containing the elements of `basis`. This is the
parent basis for a [`SubBasis`](@ref) and `basis` itself otherwise.
"""
implicit_basis(basis::AbstractBasis) = basis
implicit_basis(basis::SubBasis) = parent(basis)

function MA.promote_operation(
    ::typeof(implicit_basis),
    ::Type{B},
) where {B<:AbstractBasis}
    return B
end

function MA.promote_operation(
    ::typeof(implicit_basis),
    ::Type{<:SubBasis{T,I,K,B}},
) where {T,I,K,B}
    return B
end

# The object of the algebra does not depend on its basis
function _implicit_algebra(alg::StarAlgebra, basis)
    mstr, _ = promote_with_map(mstructure(alg), basis, identity)
    return StarAlgebra(object(alg), mstr)
end

"""
    implicit(a::Union{AlgebraElement,Term})

Return `a` as an element of the algebra whose basis is
`implicit_basis(basis(a))`.

The result of an arithmetic operation involving an element of the algebra of a
[`SubBasis`](@ref) is generally not in the span of this `SubBasis` so the
operands are converted with `implicit` first. This only depends on the type of
the basis so that the type of the result can be inferred.
"""
implicit(a::AlgebraElement) = _implicit(a, basis(a), implicit_basis(basis(a)))
_implicit(a::AlgebraElement, ::B, ::B) where {B} = a
function _implicit(a::AlgebraElement, _, basis)
    return AlgebraElement(coeffs(a, basis), _implicit_algebra(parent(a), basis))
end

function implicit(t::Term)
    return _implicit(t, basis(t.algebra), implicit_basis(basis(t.algebra)))
end
_implicit(t::Term, ::B, ::B) where {B} = t
function _implicit(t::Term, _, basis)
    return Term(
        _implicit_algebra(t.algebra, basis),
        basis[basis_element(t)],
        coefficient(t),
    )
end

function MA.promote_operation(::typeof(implicit), ::Type{A}) where {A}
    return Base.promote_op(implicit, A)
end

_has_sub_basis(::Type{A}) where {A} = Base.promote_op(basis, A) <: SubBasis

function _promote_operands(a, b)
    return promote_bases(implicit(a), implicit(b))
end

# Parents can have different types when their basis data is encoded in types.
# Reuse the same promotion as values before applying coefficient promotion.
function _promoted_types(::Type{A}, ::Type{B}) where {A,B}
    P = Base.promote_op(_promote_operands, A, B)
    return isconcretetype(P) ? fieldtypes(P) : (Any, Any)
end
function _convert_parent(::Type{A}, x) where {A<:AbstractStarAlgebra}
    target = convert(A, parent(x))
    y, _ = promote_bases(x, target)
    parent(y) == target || throw(InexactError(:convert, A, parent(x)))
    return y
end
function _zero_prototype(x::AlgebraElement, ::Type{T}) where {T}
    return MA.operate!(zero, similar(x, T))
end
_zero_prototype(x::Term, ::Type{T}) where {T} = zero(T, parent(x))

@inline function _zero_for(::Type{P}, alg, prototype) where {P<:AlgebraElement}
    isnothing(prototype) && return zero(P)
    z = _zero_prototype(prototype, eltype(P))
    return convert(P, first(promote_bases(z, alg)))
end

# Scalar representations can expose their term for parent preparation.
_algebra_operand(x) = x

_promote_parent(::Nothing, ::Nothing) = nothing
_promote_parent(::Nothing, b) = b
_promote_parent(a, ::Nothing) = a
_promote_parent(a, b) = first(promote_bases(a, b))

"""
    common_parent(operands)

Return the parent to which the parents of the algebra elements and terms of the
collection `operands` can all be promoted. If it returns `nothing`, which is the
default, they are promoted one by one with [`promote_bases`](@ref). It can be
specialized to promote all of them at once.
"""
common_parent(_) = nothing

# Keep optional parents and prototypes unboxed across input recursion.
@inline function _zero_for(
    ::Type{P},
    alg,
    prototype,
    input,
    inputs::Vararg{Any,N},
) where {P<:AlgebraElement,N}
    common = common_parent(input)
    alg = _promote_parent(alg, common)
    for operand in input
        x = _algebra_operand(operand)
        if x isa Union{AlgebraElement,Term}
            if isnothing(common)
                alg = _promote_parent(alg, parent(x))
            end
            # Prefer an element with the promoted coefficient storage, e.g.
            # dense storage when mixing sparse and dense coefficients.
            if isnothing(prototype) || (
                x isa AlgebraElement &&
                similar_type(typeof(x), eltype(P)) === P
            )
                prototype = x
            end
        end
    end
    # Recurse over the tuple so each input's iterator type stays concrete.
    return _zero_for(P, alg, prototype, inputs...)
end

function MA.zero_for(
    ::Type{P},
    inputs::Vararg{Any,N},
) where {P<:AlgebraElement,N}
    return _zero_for(P, nothing, nothing, inputs...)
end
