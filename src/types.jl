# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

abstract type AbstractStarAlgebra{O,T} end

mstructure(A::AbstractStarAlgebra) = mstructure(basis(A))

function _sanity_checks(coeffs, A::AbstractStarAlgebra)
    @assert key_type(coeffs) == key_type(basis(A))
end
function _sanity_checks(coeffs::AbstractVector, A::AbstractStarAlgebra)
    @assert Base.haslength(basis(A))
    @assert length(coeffs) == length(basis(A))
end

# concrete implementation
struct StarAlgebra{O,T,M<:MStructure{T}} <: AbstractStarAlgebra{O,T}
    object::O
    mstructure::M
end

mstructure(A::StarAlgebra) = A.mstructure
basis(A::StarAlgebra) = basis(mstructure(A))
function MA.promote_operation(::typeof(basis), ::Type{StarAlgebra{O,T,M}}) where {O,T,M}
    return MA.promote_operation(bass, M)
end
object(A::StarAlgebra) = A.object

struct AlgebraElement{A,T,V} <: MA.AbstractMutable
    coeffs::V
    parent::A
end

Base.parent(a::AlgebraElement) = a.parent
Base.eltype(::Type{A}) where {A<:AlgebraElement} = value_type(MA.promote_operation(coeffs, A))
Base.eltype(a::AlgebraElement) = eltype(typeof(a))
function MA.promote_operation(::typeof(coeffs), ::Type{AlgebraElement{A,T,V}}) where {A,T,V}
    return V
end
coeffs(a::AlgebraElement) = a.coeffs

function MA.operate!(T::typeof(canonical), a::AlgebraElement)
    return MA.operate!(T, coeffs(a))
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
Base.zero(A::AbstractStarAlgebra) = zero(Int, A)
function Base.zero(T::Type, A::AbstractStarAlgebra)
    cfs = zero_coeffs(T, basis(A))
    return AlgebraElement(cfs, A)
end

Base.one(A::AbstractStarAlgebra) = one(Int, A)
function Base.one(T::Type, A::AbstractStarAlgebra)
    i = one(object(A))
    sc = SparseCoefficients([i], [one(T)])
    # TODO
    # this is not correct, more thought is needed
    if basis(A) isa DiracBasis
        @assert i in basis(A)
        return AlgebraElement(sc, A)
    else
        return AlgebraElement(
            coeffs(sc, DiracBasis(object(A)), basis(A)),
            A,
        )
    end
end

Base.zero(a::AlgebraElement) = (b = similar(a); return MA.operate!(zero, b))
Base.one(a::AlgebraElement) = one(eltype(a), parent(a))
Base.iszero(a::AlgebraElement) = iszero(coeffs(a))

function Base.isone(a::AlgebraElement)
    c = coeffs(a)
    A = parent(a)
    cfs1 = SparseCoefficients((one(object(A)),), (1,))

    if basis(A) isa DiracBasis
        return c == cfs1
    else
        dc = coeffs(c, basis(a), DiracBasis(object(parent(a))))
        return dc == cfs1
    end
end

function (A::AbstractStarAlgebra{O,T})(elt::T) where {O,T}
    b = basis(A)
    @assert elt in b
    res = zero_coeffs(Int, b)
    res[b[elt]] = 1
    return AlgebraElement(res, A)
end

(A::AbstractStarAlgebra)(x::Number) = x * one(A)

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
