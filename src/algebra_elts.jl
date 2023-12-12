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
