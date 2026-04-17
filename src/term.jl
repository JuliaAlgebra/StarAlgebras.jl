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

Base.:^(x::Term, p::Integer) = Base.power_by_squaring(x, p)

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

function MA.mutability(::Type{Term{T,A,I}}) where {T,A,I}
    if MA.mutability(T) isa MA.IsMutable && MA.mutability(I) isa MA.IsMutable
        return MA.IsMutable()
    else
        return MA.IsNotMutable()
    end
end

# dot for Term to avoid recursive fallback in LinearAlgebra.dot
function LinearAlgebra.dot(t1::Term, t2::Term)
    return coefficient(t1) *
           coefficient(t2) *
           (basis_element(t1) * basis_element(t2))
end
function LinearAlgebra.dot(x, t::Term)
    return x * t
end
function LinearAlgebra.dot(t::Term, x)
    return star(t) * x
end

function MA.operate_to!(t::Term, ::typeof(*), t1::Term, t2::Term)
    MA.operate_to!(t.coefficient, *, coefficient(t1), coefficient(t2))
    # Index multiplication goes through the algebra's mstructure
    ms = mstructure(t.algebra)
    # For the index, we can't mutate in general, so reconstruct
    # TODO: add mutable index path for performance
    new_idx = ms(t1.index, t2.index, eltype(basis(t.algebra)))
    # new_idx is a SparseCoefficients with one entry for commutative case
    # Extract the single key
    @assert length(collect(keys(new_idx))) == 1 "Term * Term must produce a single basis element"
    t_new = Term(t.algebra, first(keys(new_idx)), t.coefficient)
    return t_new
end

function MA.operate!(::typeof(*), t1::Term, t2::Term)
    MA.operate!(*, t1.coefficient, coefficient(t2))
    ms = mstructure(t1.algebra)
    new_idx = ms(t1.index, t2.index, eltype(basis(t1.algebra)))
    @assert length(collect(keys(new_idx))) == 1
    # Can't mutate index field of immutable struct, return new term
    return Term(t1.algebra, first(keys(new_idx)), t1.coefficient)
end

function MA.operate!(::typeof(one), t::Term)
    MA.operate!(one, t.coefficient)
    # Can't mutate index to identity in immutable struct
    # Return new term with identity index
    return Term(t.algebra, t.index, t.coefficient)
end

# Term + Term → AlgebraElement
function Base.:+(t1::Term, t2::Term)
    (a1, m1), (a2, m2) = promote_bases_with_maps(t1.algebra, t2.algebra)
    idx1 = m1 === nothing ? t1.index : m1(t1.index)
    idx2 = m2 === nothing ? t2.index : m2(t2.index)
    sc = SparseCoefficients([idx1, idx2], [coefficient(t1), coefficient(t2)])
    return algebra_element(sc, a1)
end

function Base.:-(t1::Term, t2::Term)
    (a1, m1), (a2, m2) = promote_bases_with_maps(t1.algebra, t2.algebra)
    idx1 = m1 === nothing ? t1.index : m1(t1.index)
    idx2 = m2 === nothing ? t2.index : m2(t2.index)
    sc = SparseCoefficients([idx1, idx2], [coefficient(t1), -coefficient(t2)])
    return algebra_element(sc, a1)
end

function Base.:-(t::Term)
    return Term(t.algebra, t.index, -coefficient(t))
end

# Term + AlgebraElement and vice versa
function Base.:+(t::Term, a::AlgebraElement)
    return algebra_element(t) + a
end
function Base.:+(a::AlgebraElement, t::Term)
    return a + algebra_element(t)
end
function Base.:-(t::Term, a::AlgebraElement)
    return algebra_element(t) - a
end
function Base.:-(a::AlgebraElement, t::Term)
    return a - algebra_element(t)
end

# Convert Term to single-entry AlgebraElement
function algebra_element(t::Term)
    sc = SparseCoefficients((t.index,), (coefficient(t),))
    return algebra_element(sc, t.algebra)
end
