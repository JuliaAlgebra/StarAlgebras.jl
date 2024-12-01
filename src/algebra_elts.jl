Base.hash(a::AlgebraElement, h::UInt) = hash(coeffs(a), hash(parent(a), h))

function Base.:(==)(X::AlgebraElement, Y::AlgebraElement)
    parent(X) === parent(Y) || return false
    return coeffs(X) == coeffs(Y)
end

Base.copy(a::AlgebraElement) = AlgebraElement(copy(coeffs(a)), parent(a))

function Base.deepcopy_internal(a::AlgebraElement, id::IdDict)
    if !haskey(id, a)
        id[a] = AlgebraElement(Base.deepcopy_internal(coeffs(a), id), parent(a))
    end
    return id[a]
end

# call overload:
(a::AlgebraElement)(x) = coeffs(a)[basis(a)[x]]
Base.setindex!(a::AlgebraElement, v, idx) = a.coeffs[basis(a)[idx]] = v

# AlgebraElement specific functions

function supp(a::AlgebraElement)
    b = basis(a)
    return [b[i] for (i, _) in nonzero_pairs(coeffs(a))]
end

function LinearAlgebra.norm(a::AlgebraElement, p::Real)
    return LinearAlgebra.norm(coeffs(a), p)
end

function LinearAlgebra.dot(a::AlgebraElement, b::AlgebraElement)
    return LinearAlgebra.dot(coeffs(a), coeffs(b))
end

function LinearAlgebra.dot(w::AbstractVector, b::AlgebraElement)
    return LinearAlgebra.dot(w, coeffs(b))
end
function LinearAlgebra.dot(a::AlgebraElement, w::AbstractVector)
    @assert key_type(basis(parent(a))) <: Integer
    return LinearAlgebra.dot(coeffs(a), w)
end
