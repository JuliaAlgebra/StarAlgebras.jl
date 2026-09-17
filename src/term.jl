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

function Base.convert(::Type{Term{T,A,I}}, t::Term{<:Any,A,I}) where {T,A,I}
    return Term(t.algebra, t.index, convert(T, coefficient(t)))
end
Base.convert(::Type{Term{T,A,I}}, t::Term{T,A,I}) where {T,A,I} = t

Base.:^(x::Term, p::Integer) = algebra_element(x)^p

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

for op in (:+, :-, :*)
    @eval begin
        Base.$op(t1::Term, t2::Term) =
            $op(algebra_element(t1), algebra_element(t2))
        Base.$op(t::Term, a::AlgebraElement) = $op(algebra_element(t), a)
        Base.$op(a::AlgebraElement, t::Term) = $op(a, algebra_element(t))
    end
end

# Convert Term to single-entry AlgebraElement
function algebra_element(t::Term)
    c = zero_coeffs(typeof(coefficient(t)), basis(parent(t)))
    if !iszero(t)
        c[t.index] = coefficient(t)
    end
    return AlgebraElement(c, parent(t))
end

function Base.convert(
    ::Type{AlgebraElement{T,A,C}},
    t::Term{S,A,I},
) where {T,A,C,S,I}
    return convert(AlgebraElement{T,A,C}, algebra_element(t))
end
