star(x::Real) = x
star(x::Complex) = conj(x)
# star(x) = x' # ???
star(::AbstractBasis, x) = star(x)
star(basis::AbstractBasis, i::Integer) = basis[-i]

function star(basis::AbstractBasis, d::Dirac)
    return Dirac(star(basis, d.element), star(d.value))
end
function star(basis::AbstractBasis, ad::AugmentedDirac)
    return AugmentedDirac(star(basis, ad.dirac))
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
    dropzeros!(v)
    return v
end
