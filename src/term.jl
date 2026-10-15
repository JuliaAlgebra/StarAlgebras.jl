# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

"""
    struct Term{T,A<:AbstractStarAlgebra,I} <: MA.AbstractMutable
        algebra::A
        index::I
        coefficient::T
    end

A term in a star algebra: a scalar coefficient times a basis element.
The basis element is identified by its index in the algebra's basis,
and can be reconstructed via `basis(algebra)[index]`.

Having the algebra as a field allows generic arithmetic:
`+(::Term, ::Term)` returns an `AlgebraElement` by combining both terms
in a common algebra (via `promote_bases`).
`*(::Term, ::Term)` returns a `Term` for [`OrderedTermProduct`](@ref),
and an `AlgebraElement` for [`GeneralTermProduct`](@ref).
Nonnegative integer powers follow the same rule.
"""
struct Term{T,A<:AbstractStarAlgebra,I} <: MA.AbstractMutable
    algebra::A
    index::I
    coefficient::T
end

# Copy constructor (needed e.g. for array operations)
Term{T,A,I}(t::Term{T,A,I}) where {T,A,I} = t

"""
    coefficient(t::Term)

Return the coefficient of the term `t`.
"""
coefficient(t::Term) = t.coefficient

"""
    basis_element(t::Term)

Return the basis element of the term `t`, reconstructed from the algebra and index.
"""
basis_element(t::Term) = basis(t.algebra)[t.index]

Base.parent(t::Term) = t.algebra

Base.iszero(t::Term) = iszero(coefficient(t))
Base.isone(t::Term) = isone(coefficient(t)) && isone(basis_element(t))

Base.zero(t::Term) = Term(t.algebra, t.index, zero(coefficient(t)))

function Base.one(t::Term)
    A = parent(t)
    return Term(A, basis(A)[one(object(A))], one(coefficient(t)))
end

function MA.promote_operation(::typeof(one), ::Type{Term{T,A,I}}) where {T,A,I}
    return Term{MA.promote_operation(one, T),A,I}
end

function remove_leading_term(t::Term)
    return zero(t)
end
function MA.operate(::typeof(remove_leading_term), t::Term)
    return remove_leading_term(t)
end
function MA.promote_operation(
    ::typeof(remove_leading_term),
    ::Type{Term{T,A,I}},
) where {T,A,I}
    return Term{MA.promote_operation(zero, T),A,I}
end

function Base.:(==)(t1::Term, t2::Term)
    c1 = coefficient(t1)
    c2 = coefficient(t2)
    if iszero(c1)
        return iszero(c2)
    end
    c1 == c2 || return false
    if t1.algebra === t2.algebra
        return t1.index == t2.index
    else
        return basis_element(t1) == basis_element(t2)
    end
end
function Base.isequal(t1::Term, t2::Term)
    c1 = coefficient(t1)
    c2 = coefficient(t2)
    if iszero(c1)
        return iszero(c2)
    end
    isequal(c1, c2) || return false
    if t1.algebra === t2.algebra
        return isequal(t1.index, t2.index)
    else
        return isequal(basis_element(t1), basis_element(t2))
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

function Base.convert(::Type{Term{T,A,I}}, t::Term{S,B,J}) where {T,A,I,S,B,J}
    y = A === B ? t : _convert_parent(A, t)
    return Term{T,A,I}(
        parent(y),
        convert(I, y.index),
        convert(T, coefficient(y)),
    )
end
Base.convert(::Type{Term{T,A,I}}, t::Term{T,A,I}) where {T,A,I} = t

function Base.promote_rule(
    ::Type{Term{T,A,I}},
    ::Type{Term{S,B,J}},
) where {T,S,A,I,B,J}
    if A !== B
        return promote_type(_promoted_types(Term{T,A,I}, Term{S,B,J})...)
    end
    return Term{promote_type(T, S),A,promote_type(I, J)}
end

function Base.promote_rule(
    ::Type{P},
    ::Type{Term{S,B,I}},
) where {T,S,A,C,B,I,P<:AlgebraElement{T,A,C}}
    if A !== B
        return promote_type(_promoted_types(P, Term{S,B,I})...)
    end
    return similar_type(P, promote_type(T, S))
end

# `power_by_squaring` promotes to the type of `x * x`, so the product style
# also determines whether a power is a Term or an AlgebraElement.
Base.:^(x::Term, p::Integer) = Base.power_by_squaring(x, p)

function MA.promote_operation(
    ::typeof(^),
    ::Type{T},
    ::Type{<:Integer},
) where {T<:Term}
    return MA.promote_operation(*, T, T)
end

Base.ndims(::Union{Type{<:Term},Term}) = 0
Base.broadcastable(t::Term) = Ref(t)

Base.copy(t::Term) = MA.mutable_copy(t)
function MA.mutable_copy(t::Term)
    return Term(
        t.algebra,
        MA.copy_if_mutable(t.index),
        MA.copy_if_mutable(coefficient(t)),
    )
end

function LinearAlgebra.dot(t1::Term, t2::Term)
    return LinearAlgebra.dot(algebra_element(t1), algebra_element(t2))
end

function Base.:-(t::Term)
    return Term(t.algebra, t.index, -coefficient(t))
end

Base.:+(x::Union{AlgebraElement,Term}) = x
Base.:*(x::Union{AlgebraElement,Term}) = x

function Base.:*(a::Union{T,Number}, t::Term{T}) where {T}
    return Term(parent(t), MA.copy_if_mutable(t.index), a * coefficient(t))
end

function Base.:*(t::Term{T}, a::Union{T,Number}) where {T}
    return Term(parent(t), MA.copy_if_mutable(t.index), coefficient(t) * a)
end

function MA.promote_operation(
    ::typeof(*),
    ::Type{S},
    ::Type{Term{T,A,I}},
) where {T,A,I,S<:Union{T,Number}}
    return Term{MA.promote_operation(*, S, T),A,I}
end

function MA.promote_operation(
    ::typeof(*),
    ::Type{Term{T,A,I}},
    ::Type{S},
) where {T,A,I,S<:Union{T,Number}}
    return Term{MA.promote_operation(*, T, S),A,I}
end

Base.:*(t::Term, a::AlgebraElement) = algebra_element(t) * a
Base.:*(a::AlgebraElement, t::Term) = a * algebra_element(t)

_addition_coefficients(a::AlgebraElement) = coeffs(a)
function _addition_coefficients(t::Term)
    b = basis(parent(t))
    if b isa ImplicitBasis
        return SparseCoefficients((t.index,), (coefficient(t),), comparable(b))
    end
    return coeffs(algebra_element(t))
end

function _term_add(op::F, a, b) where {F}
    a, b = promote_bases(a, b)
    P = MA.promote_operation(op, typeof(a), typeof(b))
    output = MA.zero_for(P, (a,), (b,))
    # Use read-only term coefficients in the existing binary kernel so that
    # coefficient arithmetic happens before conversion to the result type.
    MA.operate_to!(
        coeffs(output),
        op,
        _addition_coefficients(a),
        _addition_coefficients(b),
    )
    return output
end

for (L, R) in ((Term, Term), (Term, AlgebraElement), (AlgebraElement, Term)),
    op in (:+, :-)

    @eval Base.$op(a::$L, b::$R) = _term_add($op, a, b)
end

function _constant_term(p, scalar)
    alg = parent(p)
    return Term(alg, basis(alg)[one(object(alg))], scalar)
end

for op in (:+, :-)
    @eval begin
        Base.$op(t::Term{T}, a::Union{T,Number}) where {T} =
            $op(t, _constant_term(t, a))
        Base.$op(a::Union{T,Number}, t::Term{T}) where {T} =
            $op(_constant_term(t, a), t)
    end
end

function Base.:*(a::Term, b::Term)
    a, b = promote_bases(a, b)
    ms = mstructure(parent(a))
    style = term_product_style(ms, comparable(basis(parent(a))))
    return _term_product(style, ms, a, b)
end

function _term_product(::OrderedTermProduct, ms, a::Term, b::Term)
    key, weight = only(nonzero_pairs(_basis_product(ms, a.index, b.index)))
    return Term(
        parent(a),
        MA.copy_if_mutable(key),
        weight * (coefficient(a) * coefficient(b)),
    )
end

function _term_product(::GeneralTermProduct, ms, a::Term, b::Term)
    return algebra_element(a) * algebra_element(b)
end

function MA.promote_operation(
    ::typeof(*),
    ::Type{A},
    ::Type{B},
) where {A<:Term,B<:Term}
    return Base.promote_op(*, A, B)
end

for (L, R, ops) in (
        (Term, Term, (+, -)),
        (Term, AlgebraElement, (+, -, *)),
        (AlgebraElement, Term, (+, -, *)),
    ),
    op in ops

    @eval function MA.promote_operation(
        ::typeof($op),
        ::Type{A},
        ::Type{B},
    ) where {A<:$L,B<:$R}
        return MA.promote_operation(
            $op,
            _algebra_element_type(A),
            _algebra_element_type(B),
        )
    end
end

_algebra_element_type(::Type{A}) where {A<:AlgebraElement} = A
function _algebra_element_type(::Type{T}) where {T<:Term}
    return Base.promote_op(algebra_element, T)
end

# Convert Term to single-entry AlgebraElement
function algebra_element(t::Term{T}) where {T}
    c = zero_coeffs(T, basis(parent(t)))
    if !iszero(t)
        c[t.index] = coefficient(t)
    end
    return AlgebraElement(c, parent(t))
end

function Base.convert(
    ::Type{AlgebraElement{T,A,C}},
    t::Term{S,B,I},
) where {T,A,C,S,B,I}
    y = A === B ? t : _convert_parent(A, t)
    return convert(AlgebraElement{T,A,C}, algebra_element(y))
end
