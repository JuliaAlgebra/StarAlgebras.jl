# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

Base.adjoint(a::AlgebraElement) = star(a)
star(x::Any) = x'

function star(X::AlgebraElement)
    res = star(mstructure(X), coeffs(X))
    return AlgebraElement(res, parent(X))
end

star(::MultiplicativeStructure, x) = star(x)

function star(mstr::MultiplicativeStructure, d::SparseCoefficients)
    k = star.(Ref(mstr), keys(d))
    v = star.(values(d))
    return SparseCoefficients(k, v)
end

function star(mstr::MTable{T}, el::T) where T
    return basis(mstr)[mstr.starof[basis(mstr)[el]]]
end

function star(mstr::MTable, coeffs::SparseVector)
    nzidcs = mstr.starof[SparseArrays.nonzeroinds(coeffs)]
    nzvals = star.(SparseArrays.nonzeros(coeffs))

    return SparseVector(length(coeffs), nzidcs, nzvals)
end
