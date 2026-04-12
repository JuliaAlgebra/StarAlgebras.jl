# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

"""
    struct Term{T,B}
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

Base.zero(t::Term) = Term(zero(coefficient(t)), basis_element(t))

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
