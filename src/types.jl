# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

abstract type AbstractStarAlgebra{O,T} end

function _sanity_checks(coeffs, A::AbstractStarAlgebra)
    @assert key_type(coeffs) == key_type(basis(A))
end
function _sanity_checks(coeffs::AbstractVector, A::AbstractStarAlgebra)
    @assert Base.haslength(basis(A))
    @assert length(coeffs) == length(basis(A))
end

"""
    struct StarAlgebra{O,T,M<:MultiplicativeStructure{T}} <: AbstractStarAlgebra{O,T}
        object::O
        mstructure::M
    end

Star algebra implementation with an `object` that should implement `one(::O)` and
a [`MultiplicativeStructure`](@ref) `mstructure`.
"""
struct StarAlgebra{O,T,M<:MultiplicativeStructure{T}} <: AbstractStarAlgebra{O,T}
    object::O
    mstructure::M
end

StarAlgebra(object, basis::AbstractBasis) = StarAlgebra(object, DiracMStructure(basis, *))

mstructure(A::StarAlgebra) = A.mstructure
basis(A::StarAlgebra) = basis(mstructure(A))
function MA.promote_operation(::typeof(basis), ::Type{StarAlgebra{O,T,M}}) where {O,T,M}
    return MA.promote_operation(basis, M)
end
object(A::StarAlgebra) = A.object
Base.isempty(A::StarAlgebra) = isempty(basis(A))

struct AlgebraElement{A,T,V} <: MA.AbstractMutable
    coeffs::V
    parent::A
end

Base.parent(a::AlgebraElement) = a.parent
Base.in(x::AlgebraElement, A::AbstractStarAlgebra) = parent(x) == A

mstructure(a::AlgebraElement) = mstructure(parent(a))
Base.eltype(::Type{A}) where {A<:AlgebraElement} = value_type(MA.promote_operation(coeffs, A))
Base.eltype(a::AlgebraElement) = eltype(typeof(a))
function MA.promote_operation(::typeof(coeffs), ::Type{AlgebraElement{A,T,V}}) where {A,T,V}
    return V
end
coeffs(a::AlgebraElement) = a.coeffs

comparable(::Type) = isless
comparable(object) = comparable(eltype(object))
key_isless(basis::AbstractBasis) = comparable(object(basis))
key_isless(mstr::MultiplicativeStructure) = key_isless(basis(mstr))

function MA.operate!(T::typeof(canonical), a::AlgebraElement)
    return MA.operate!(T, coeffs(a), key_isless(basis(a)))
end

function coeffs(x::AlgebraElement, b::AbstractBasis)
    return coeffs(coeffs(x), basis(x), b)
end
function adjoint_coeffs(a::AlgebraElement, target::AbstractBasis)
    return adjoint_coeffs(coeffs(a), target, basis(a))
end
function MA.promote_operation(::typeof(basis), ::Type{<:AlgebraElement{A}}) where {A}
    return MA.promote_operation(basis, A)
end
basis(a::AlgebraElement) = basis(parent(a))

function AlgebraElement(coeffs, A::AbstractStarAlgebra)
    _sanity_checks(coeffs, A)
    return AlgebraElement{typeof(A),value_type(coeffs),typeof(coeffs)}(coeffs, A)
end

function AlgebraElement(
    coeffs::SparseCoefficients{T},
    A::AbstractStarAlgebra{O,T},
) where {O,T}
    return AlgebraElement{typeof(A),value_type(coeffs),typeof(coeffs)}(coeffs, A)
end

### constructing elements
function __coerce(A::AbstractStarAlgebra, (x,v)::Pair{K, V}) where {K,V}
    if iszero(v)
        return AlgebraElement(zero_coeffs(V, basis(A)), A)
    elseif x in basis(A)
        cfs = zero_coeffs(V, basis(A))
        setindex!(cfs, v, basis(A)[x]; lt = key_isless(basis(A)))
        return AlgebraElement(cfs, A)
    # elseif x in object(A)
    #     sc = SparseCoefficients([x], [v])
    #     return AlgebraElement(
    #         coeffs(sc, DiracBasis(object(A)), basis(A)),
    #         A,
    #     )
    else
        throw(ArgumentError("cannot coerce $x to $A"))
    end
end

Base.zero(A::AbstractStarAlgebra) = zero(Int, A)
Base.zero(T::Type, A::AbstractStarAlgebra) =
    __coerce(A, (one(object(A)) => zero(T)))
Base.zero(a::AlgebraElement) = (b = similar(a); return MA.operate!(zero, b))
Base.iszero(a::AlgebraElement) = iszero(coeffs(a))

Base.one(A::AbstractStarAlgebra) = one(Int, A)
Base.one(T::Type, A::AbstractStarAlgebra) =
    __coerce(A, (one(object(A)) => one(T)))
Base.one(a::AlgebraElement) = one(eltype(a), parent(a))

function Base.isone(a::AlgebraElement)
    A = parent(a)
    id = one(object(A))
    if id in basis(A)
        for (i, v) in nonzero_pairs(coeffs(a))
            isone(v) || return false
            id == basis(A)[i] || return false
        end
        return true
    else
        throw(ArgumentError("basis of $A does not contain $id; `one` and `isone` are unsupported"))
        # return a == one(a)
    end
end

(A::AbstractStarAlgebra{O,T})(elt::T) where {O,T} = __coerce(A, (elt => 1))
(A::AbstractStarAlgebra)(x::Number) = __coerce(A, (one(object(A)) => x))

function similar_type(::Type{AlgebraElement{A,T,V}}, ::Type{C}) where {A,T,V,C}
    return AlgebraElement{A,C,similar_type(V, C)}
end

function Base.similar(X::AlgebraElement, T = eltype(X))
    return AlgebraElement(similar(coeffs(X), T), parent(X))
end

function AlgebraElement{T}(X::AlgebraElement) where {T}
    v = coeffs(X)
    w = similar(v, T)
    MA.operate!(zero, w)
    for (k, v) in nonzero_pairs(v)
        w[k] = v
    end
    return AlgebraElement(w, parent(X))
end
