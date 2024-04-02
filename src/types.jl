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
struct StarAlgebra{O,T,B<:AbstractBasis{T}} <:
       AbstractStarAlgebra{O,T}
    object::O
    basis::B

    function StarAlgebra(obj, basis::AbstractBasis)
        O = typeof(obj)
        T = eltype(basis)
        B = typeof(basis)

        return new{O,T,B}(obj, basis)
    end
end

# MTable:
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

basis(A::StarAlgebra) = A.basis
object(A::StarAlgebra) = A.object

struct AlgebraElement{A,T,V} <: MA.AbstractMutable
    coeffs::V
    parent::A
end

Base.parent(a::AlgebraElement) = a.parent
Base.eltype(a::AlgebraElement) = valtype(coeffs(a))
coeffs(a::AlgebraElement) = a.coeffs
function coeffs(x::AlgebraElement, b::AbstractBasis)
    return coeffs(coeffs(x), basis(parent(x)), b)
end

function AlgebraElement(coeffs, A::AbstractStarAlgebra)
    _sanity_checks(coeffs, A)
    return AlgebraElement{typeof(A),valtype(coeffs),typeof(coeffs)}(coeffs, A)
end

function AlgebraElement(
    coeffs::SparseCoefficients{T},
    A::AbstractStarAlgebra{O,T}
) where {O,T}
    return AlgebraElement{typeof(A),valtype(coeffs),typeof(coeffs)}(coeffs, A)
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
            coeffs(sc, DiracBasis{UInt}(object(A)), basis(A)),
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
        dc = coeffs(c, basis(parent(a)), DiracBasis{UInt}(object(parent(a))))
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

Base.similar(X::AlgebraElement, T=eltype(X)) = AlgebraElement(similar(coeffs(X), T), parent(X))

function AlgebraElement{T}(X::AlgebraElement) where {T}
    v = coeffs(X)
    w = similar(v, T)
    MA.operate!(zero, w)
    for (k, v) in nonzero_pairs(v)
        w[k] = v
    end
    return AlgebraElement(w, parent(X))
end
