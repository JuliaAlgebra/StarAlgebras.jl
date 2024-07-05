Base.adjoint(a::AlgebraElement) = star(a)
star(x::Any) = x'

function star(X::AlgebraElement)
    res = star(basis(X), coeffs(X))
    return AlgebraElement(res, parent(X))
end

star(::AbstractBasis, x) = star.(x)

function star(basis::AbstractBasis, d::SparseCoefficients)
    k = star.(Ref(basis), keys(d))
    v = star.(values(d))
    return SparseCoefficients(k, v)
end

function star(basis::FixedBasis, coeffs::SparseVector)
    nzidcs = mstructure(basis).starof[SparseArrays.nonzeroinds(coeffs)]
    nzvals = star.(SparseArrays.nonzeros(coeffs))

    v = SparseVector(length(coeffs), nzidcs, nzvals)
    return v
end
