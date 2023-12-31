Base.hash(a::AlgebraElement, h::UInt) = hash(coeffs(a), hash(parent(a), h))

function Base.:(==)(X::AlgebraElement, Y::AlgebraElement)
    parent(X) === parent(Y) || return false
    return coeffs(X) == coeffs(Y)
end

Base.getindex(a::AlgebraElement, x) = coeffs(a)[x]
Base.setindex!(a::AlgebraElement, v, idx) = a.coeffs[idx] = v
# call overload:
(a::AlgebraElement)(x) = a[x]

# AlgebraElement specific functions

supp(a::AlgebraElement) = (b = basis(parent(a)); [b[i] for i in keys(coeffs(a))])

function star(X::AlgebraElement)
    res = star(basis(parent(X)), coeffs(X))
    return AlgebraElement(res, parent(X))
end

Base.adjoint(a::AlgebraElement) = star(a)

LinearAlgebra.norm(a::AlgebraElement, p::Real) =
    LinearAlgebra.norm(coeffs(a), p)
aug(a::AlgebraElement) = sum(coeffs(a))

LinearAlgebra.dot(a::AlgebraElement, v::AbstractVector) =
    LinearAlgebra.dot(StarAlgebras.coeffs(a), v)
LinearAlgebra.dot(v::AbstractVector, a::AlgebraElement) = LinearAlgebra.dot(a, v)

Base.copy(a::AlgebraElement) = AlgebraElement(copy(coeffs(a)), parent(a))

function Base.deepcopy_internal(a::AlgebraElement, IdDict::IdDict)
    if !haskey(IdDict, a)
        IdDict[a] = AlgebraElement(Base.deepcopy_internal(coeffs(a), IdDict), parent(a))
    end
    return IdDict[a]
end
