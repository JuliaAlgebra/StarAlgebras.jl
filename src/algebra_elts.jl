Base.hash(a::AlgebraElement, h::UInt) = hash(coeffs(a), hash(parent(a), h))

function Base.:(==)(X::AlgebraElement, Y::AlgebraElement)
    parent(X) === parent(Y) || return false
    return coeffs(X) == coeffs(Y)
end

Base.getindex(a::AlgebraElement, i::Integer) = coeffs(a)[i]
Base.getindex(a::AlgebraElement{<:StarAlgebra{O,T}}, x::T) where {O,T} =
    (b = basis(parent(a)); a[b[x]])

# call overload:
(a::AlgebraElement)(x) = a[x]

Base.setindex!(a::AlgebraElement, v, i::Integer) = a.coeffs[i] = v

function Base.setindex!(
    a::AlgebraElement{<:StarAlgebra{O,T}},
    v,
    t::T,
) where {O,T}
    b = basis(parent(a))
    return a[b[t]] = v
end

# AlgebraElement specific functions

supp_ind(a::AlgebraElement) = findall(!iszero, coeffs(a))
supp_ind(a::AlgebraElement{A,T,<:SparseVector}) where {A,T} =
    (dropzeros!(coeffs(a)); SparseArrays.nonzeroinds(coeffs(a)))
supp(a::AlgebraElement) = (b = basis(parent(a)); [b[i] for i in supp_ind(a)])

function star(A::StarAlgebra, i::Integer)
    @assert i > 0
    if i < max(size(A.mstructure)...)
        return _get(A.mstructure, -signed(i))
    else
        b = basis(A)
        return b[star(b[i])]
    end
end

function star(X::AlgebraElement)
    A = parent(X)
    b = basis(A)
    supp_X = supp_ind(X)
    idcs = similar(supp_X)
    vals = similar(idcs, eltype(X))
    for (i, idx) in enumerate(supp_X)
        idcs[i] = star(parent(X), idx)
        vals[i] = X[idx]
    end
    return AlgebraElement(sparsevec(idcs, vals, length(b)), A)
end

Base.adjoint(a::AlgebraElement) = star(a)

LinearAlgebra.norm(a::AlgebraElement, p::Real) =
    LinearAlgebra.norm(coeffs(a), p)
aug(a::AlgebraElement) = sum(coeffs(a))

LinearAlgebra.dot(a::AlgebraElement, v::AbstractVector) =
    LinearAlgebra.dot(StarAlgebras.coeffs(a), v)
LinearAlgebra.dot(v::AbstractVector, a::AlgebraElement) = LinearAlgebra.dot(a, v)

Base.copy(a::AlgebraElement) = AlgebraElement(copy(coeffs(a)), parent(a))
function Base.deepcopy_internal(a::AlgebraElement, stackdict::IdDict)
    haskey(stackdict, a) && return stackdict[a]
    return AlgebraElement(Base.deepcopy_internal(coeffs(a), stackdict), parent(a))
end
