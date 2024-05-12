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
(a::AlgebraElement)(x) = coeffs(a)[basis(parent(a))[x]]
Base.setindex!(a::AlgebraElement, v, idx) = a.coeffs[basis(parent(a))[idx]] = v

# AlgebraElement specific functions

function supp(a::AlgebraElement)
    b = basis(parent(a))
    return [b[i] for (i, _) in nonzero_pairs(coeffs(a))]
end

function LinearAlgebra.norm(a::AlgebraElement, p::Real)
    return LinearAlgebra.norm(coeffs(a), p)
end
function LinearAlgebra.dot(a::AlgebraElement, v::AbstractVector)
    return LinearAlgebra.dot(coeffs(a), v)
end
function LinearAlgebra.dot(v::AbstractVector, a::AlgebraElement)
    return LinearAlgebra.dot(a, v)
end
