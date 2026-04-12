# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

"""
    struct Term{T,B} <: MA.AbstractMutable
        coefficient::T
        basis_element::B
    end

A representation of the multiplication between a `coefficient` and a `basis_element`.
This is the generic analog of a single term in an algebra: a scalar times a basis element.

The `coefficient` does not need to be a `Number`. For instance, in multivariate
polynomial GCD computations, the coefficient can itself be a polynomial.
"""
struct Term{T,B} <: MA.AbstractMutable
    coefficient::T
    basis_element::B
end

# Copy constructor (needed e.g. for array operations)
Term{T,B}(t::Term{T,B}) where {T,B} = t

"""
    coefficient(t::Term)

Return the coefficient of the term `t`.
"""
coefficient(t::Term) = t.coefficient

"""
    basis_element(t::Term)

Return the basis element of the term `t`.
"""
basis_element(t::Term) = t.basis_element

Base.iszero(t::Term) = iszero(coefficient(t))
Base.isone(t::Term) = isone(coefficient(t)) && isone(basis_element(t))

Base.zero(t::Term) = Term(zero(coefficient(t)), basis_element(t))
Base.one(t::Term) = Term(one(coefficient(t)), one(basis_element(t)))
Base.one(::Type{Term{T,B}}) where {T,B} = Term(one(T), one(B))

function Base.:(==)(t1::Term, t2::Term)
    c1 = coefficient(t1)
    c2 = coefficient(t2)
    if iszero(c1)
        return iszero(c2)
    else
        return c1 == c2 && basis_element(t1) == basis_element(t2)
    end
end
function Base.isequal(t1::Term, t2::Term)
    c1 = coefficient(t1)
    c2 = coefficient(t2)
    if iszero(c1)
        return iszero(c2)
    else
        return isequal(c1, c2) && isequal(basis_element(t1), basis_element(t2))
    end
end

function Base.hash(t::Term, u::UInt)
    if iszero(t)
        return hash(0, u)
    elseif isone(coefficient(t))
        return hash(basis_element(t), u)
    else
        return hash(basis_element(t), hash(coefficient(t), u))
    end
end

function Base.convert(::Type{Term{T,B}}, t::Term) where {T,B}
    return Term{T,B}(convert(T, coefficient(t)), convert(B, basis_element(t)))
end
Base.convert(::Type{Term{T,B}}, t::Term{T,B}) where {T,B} = t

function Base.promote_rule(
    ::Type{Term{S,B1}},
    ::Type{Term{T,B2}},
) where {S,T,B1,B2}
    return Term{promote_type(S, T),promote_type(B1, B2)}
end

Base.:^(x::Term, p::Integer) = Base.power_by_squaring(x, p)

Base.ndims(::Union{Type{<:Term},Term}) = 0
Base.broadcastable(t::Term) = Ref(t)

# `Base.power_by_squaring` calls `Base.copy` and we want
# `t^1` to be a mutable copy of `t` so `copy` needs to be
# the same as `mutable_copy`.
Base.copy(t::Term) = MA.mutable_copy(t)
function MA.mutable_copy(t::Term)
    return Term(
        MA.copy_if_mutable(coefficient(t)),
        MA.copy_if_mutable(basis_element(t)),
    )
end

function MA.mutability(::Type{Term{T,B}}) where {T,B}
    if MA.mutability(T) isa MA.IsMutable && MA.mutability(B) isa MA.IsMutable
        return MA.IsMutable()
    else
        return MA.IsNotMutable()
    end
end

# dot for Term to avoid recursive fallback in LinearAlgebra.dot
function LinearAlgebra.dot(t1::Term, t2::Term)
    return coefficient(t1) * coefficient(t2) * (basis_element(t1) * basis_element(t2))
end
function LinearAlgebra.dot(x, t::Term)
    return x * t
end
function LinearAlgebra.dot(t::Term, x)
    return star(t) * x
end

function MA.operate_to!(
    t::Term,
    ::typeof(*),
    t1::Term,
    t2::Term,
)
    MA.operate_to!(t.coefficient, *, coefficient(t1), coefficient(t2))
    MA.operate_to!(t.basis_element, *, basis_element(t1), basis_element(t2))
    return t
end

function MA.operate!(::typeof(*), t1::Term, t2::Term)
    MA.operate!(*, t1.coefficient, coefficient(t2))
    MA.operate!(*, t1.basis_element, basis_element(t2))
    return t1
end

function MA.operate!(::typeof(one), t::Term)
    MA.operate!(one, t.coefficient)
    MA.operate!(one, t.basis_element)
    return t
end
