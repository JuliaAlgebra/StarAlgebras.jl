abstract type AbstractStarAlgebra{O,T} end

struct StarAlgebra{O,T,M<:MultiplicativeStructure,B<:AbstractBasis{T}} <:
       AbstractStarAlgebra{O,T}
    object::O
    mstructure::M
    basis::B

    function StarAlgebra(obj, basis::AbstractBasis, mstr::MultiplicativeStructure)
        O = typeof(obj)
        T = eltype(basis)
        M = typeof(mstr)
        B = typeof(basis)

        return new{O,T,M,B}(obj, mstr, basis)
    end

    function StarAlgebra(obj, mstr::MultiplicativeStructure)
        O = typeof(obj)
        T = eltype(obj)
        M = typeof(mstr)
        B = Basis{T,eltype(mstr)}

        return new{O,T,M,B}(obj, mstr)
    end
end

# TrivialMStructure:
function StarAlgebra(obj, basis::AbstractBasis)
    mstr = TrivialMStructure(basis)
    return StarAlgebra(obj, basis, mstr)
end

# CachedMStructure:
function StarAlgebra(
    obj,
    basis::AbstractBasis,
    cache_size::Tuple{<:Integer,Integer};
    precompute=false
)
    mstr = MTable(basis, size=cache_size)
    precompute && complete!(mstr)
    return StarAlgebra(obj, basis, mstr)
end

hasbasis(A::StarAlgebra) = isdefined(A, :basis)

basis(A::StarAlgebra) = A.basis
object(A::StarAlgebra) = A.object
# Base.eltype(A::StarAlgebra{O,B}) where {O,B} = eltype(B)

struct AlgebraElement{A,T,V<:AbstractVector{T}}
    coeffs::V
    parent::A

    function AlgebraElement(coeffs::AbstractVector, A::AbstractStarAlgebra)
        if hasbasis(A)
            @assert length(coeffs) == length(basis(A))
        end
        return new{typeof(A),eltype(coeffs),typeof(coeffs)}(coeffs, A)
    end
end

coeffs(a::AlgebraElement) = a.coeffs
Base.parent(a::AlgebraElement) = a.parent
Base.eltype(a::AlgebraElement) = eltype(coeffs(a))

### constructing elements
Base.zero(A::AbstractStarAlgebra) = zero(Int, A)
function Base.zero(T::Type, A::AbstractStarAlgebra)
    if hasbasis(A)
        I = SparseArrays.indtype(basis(A))
        return AlgebraElement(sparsevec(I[], T[], length(basis(A))), A)
    end
    throw(
        "Algebra without basis; use the `AlgebraElement` constructor directly.",
    )
end

Base.one(A::AbstractStarAlgebra) = one(Int, A)
function Base.one(T::Type, A::AbstractStarAlgebra)
    if hasbasis(A)
        b = basis(A)
        i = b[one(object(A))]
        return AlgebraElement(sparsevec([i], [one(T)], length(b)), A)
    end
    throw(
        "Algebra without basis; use the `AlgebraElement` constructor directly.",
    )
end

Base.zero(a::AlgebraElement) = (b = similar(a); return zero!(b))
Base.one(a::AlgebraElement) = one(parent(a))
Base.iszero(a::AlgebraElement) = iszero(coeffs(a))

function Base.isone(a::AlgebraElement)
    b = basis(parent(a))
    k = findfirst(!iszero, coeffs(a))
    k === nothing && return false
    isone(a[k]) || return false
    isone(b[k]) || return false
    return isnothing(findnext(!iszero, coeffs(a), k + 1))
end

function (A::AbstractStarAlgebra{O,T})(elt::T) where {O,T}
    if hasbasis(A)
        b = basis(A)
        i = b[elt]
        return AlgebraElement(sparsevec([i], [1], length(b)), A)
    else
        throw("Algebra without basis: cannot coerce $elt")
    end
end

function (A::AbstractStarAlgebra)(x::Number)
    if hasbasis(A)
        b = basis(A)
        i = b[one(object(A))]
        return AlgebraElement(sparsevec([i], [x], length(b)), A)
    else
        throw("Algebra without basis: cannot coerce $x")
    end
end

Base.similar(X::AlgebraElement, T=eltype(X)) = AlgebraElement(similar(coeffs(X), T), parent(X))

function AlgebraElement{T}(X::AlgebraElement) where {T}
    v = coeffs(X)
    w = similar(v, T)
    w .= v
    return AlgebraElement(w, parent(X))
end
