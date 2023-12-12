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
function Base.deepcopy_internal(a::AlgebraElement, stackdict::IdDict)
    haskey(stackdict, a) && return stackdict[a]
    return AlgebraElement(Base.deepcopy_internal(coeffs(a), stackdict), parent(a))

function fmac!(
    res::AlgebraElement{A,T,<:SparseCoefficients},
    X::AlgebraElement{A,T,<:SparseCoefficients},
    Y::AlgebraElement{A,T,<:SparseCoefficients},
) where {A,T}
    for (kx, vx) in pairs(coeffs(X))
        for (ky, vy) in pairs(coeffs(Y))
            res = MA.operate!!(MA.add_mul, res, vx * vy, basis(parent(X))[kx], basis(parent(Y))[ky])
        end
    end
    return res
end

function MA.operate!!(
    ::typeof(MA.add_mul),
    result::AlgebraElement{A,V,<:SparseCoefficients},
    α::V,
    a::DiracDelta,
    b::DiracDelta,
) where {A,V}
    (ka,), (va,) = keys(a), values(a)
    (kb,), (vb,) = keys(b), values(b)
    xy = ka * kb
    k = searchsortedfirst(coeffs(result).basis_elements, xy)
    if k in eachindex(coeffs(result).basis_elements) &&
       coeffs(result).basis_elements[k] == xy
        coeffs(result).values[k] += α * va * vb
    else
        insert!(coeffs(result).basis_elements, k, xy)
        insert!(coeffs(result).values, k, α * va * vb)
    end
    return result
end
