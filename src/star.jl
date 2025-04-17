# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

Base.adjoint(a::AlgebraElement) = star(a)
star(x::Any) = x'

function star(X::AlgebraElement)
    res = star(mstructure(X), coeffs(X))
    return AlgebraElement(res, parent(X))
end

star(mstr::MultiplicativeStructure, x) = star(basis(mstr), x)

star(::AbstractBasis, x) = star(x)

function star(b::AbstractBasis, d::SparseCoefficients)
    k = star.(Ref(b), keys(d))
    v = star.(values(d))
    return SparseCoefficients(k, v)
end

function star(b::FixedBasis, i::Integer)
    return b.starof[i]
end

function star(b::FixedBasis{T}, el::T) where T
    return b[star(b, b[el])]
end

function star(b::FixedBasis, coeffs::SparseVector)
    nzidcs = b.starof[SparseArrays.nonzeroinds(coeffs)]
    nzvals = star.(SparseArrays.nonzeros(coeffs))

    return SparseVector(length(coeffs), nzidcs, nzvals)
end
