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
        T = Symbol
        M = typeof(mstr)
        B = Basis{T,Int}

        return new{O,T,M,B}(obj, mstr)
    end
end

# TrivialMStructure:
StarAlgebra(obj, basis::AbstractBasis) = StarAlgebra{false}(obj, basis)

function StarAlgebra{Tw}(obj, basis::AbstractBasis) where {Tw}
    mstr = TrivialMStructure{Tw}(basis)
    return StarAlgebra(obj, basis, mstr)
end

# CachedMStructure:
StarAlgebra(obj, basis::AbstractBasis, cache_size::Tuple{<:Integer,Integer}; precompute=false) =
    StarAlgebra{false}(obj, basis, cache_size, precompute=precompute)

function StarAlgebra{Tw}(
    obj,
    basis::AbstractBasis,
    cache_size::Tuple{<:Integer,Integer};
    precompute=false
) where {Tw}
    mstr = CachedMTable{Tw}(basis, table_size = cache_size)
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

function Base.zero(A::AbstractStarAlgebra, T = Int)
    if hasbasis(A)
        I = SparseArrays.indtype(basis(A))
        return AlgebraElement(sparsevec(I[], T[], length(basis(A))), A)
    end
    throw(
        "Algebra without basis; to construct zero use the `AlgebraElement` constructor directly.",
    )
end

function Base.one(A::AbstractStarAlgebra, T = Int)
    hasbasis(A) && return A(one(object(A)), T)
    throw(
        "Algebra without basis; to construct one use the `AlgebraElement` constructor directly.",
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
    return isone(b[k])
end
let SA = @static VERSION < v"1.3.0" ? :StarAlgebra : :AbstractStarAlgebra
    @eval begin
        function (A::$SA{O,T})(elt::T, S = Int) where {O,T}
            if hasbasis(A)
                b = basis(A)
                i = b[elt]
                return AlgebraElement(sparsevec([i], [one(S)], length(b)), A)
            else
                throw("Algebra without basis: cannot coerce $elt.")
            end
        end

        function (A::$SA)(x::Number)
            g = one(object(A))
            res = A(g, typeof(x))
            res[g] *= x
            return res
        end
    end
end

Base.similar(X::AlgebraElement, ::Type{T} = eltype(X)) where {T} =
    AlgebraElement(similar(coeffs(X), T), parent(X))

function AlgebraElement{T}(X::AlgebraElement) where T
    v = coeffs(X)
    w = similar(v, T)
    w .= v
    return AlgebraElement(w, parent(X))
end
