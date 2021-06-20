function Base.hash(a::AlgebraElement, h::UInt)
    return hash(coeffs(a), hash(parent(a), hash(typeof(a), h)))
end
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

function Base.setindex!(a::AlgebraElement{<:StarAlgebra{O,T}}, v, t::T) where {O,T}
    b = basis(parent(a))
    return a[b[t]] = v
end


# AlgebraElement specific functions

supp_ind(a::AlgebraElement) = findall(!iszero, coeffs(a))
supp_ind(a::AlgebraElement{A,T,<:SparseVector}) where {A,T} = coeffs(a).nzind
supp(a::AlgebraElement) = (b = basis(parent(a)); [b[i] for i in supp_ind(a)])

function star(X::AlgebraElement)
    A = parent(X)
    b = basis(A)
    supp_X = supp_ind(X)
    idcs = similar(supp_X)
    vals = similar(idcs, eltype(X))
    for (i, idx) in enumerate(supp_X)
        idcs[i] = b[star(b[idx])]
        vals[i] = X[idx]
    end
    return AlgebraElement(sparsevec(idcs, vals, length(b)), A)
end

LinearAlgebra.norm(a::AlgebraElement, p) = LinearAlgebra.norm(coeffs(a), p)
aug(a::AlgebraElement) = sum(coeffs(a))
