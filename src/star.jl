star(x::Any) = x'
function star(X::AlgebraElement)
    res = star(basis(parent(X)), coeffs(X))
    return AlgebraElement(res, parent(X))
end

star(::AbstractBasis, x) = star(x)
star(basis::ImplicitBasis, i::Integer) = basis[-i]

function star(basis::AbstractBasis, ad::Augmented)
    return Augmented(star(basis, ad.elt))
end

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
