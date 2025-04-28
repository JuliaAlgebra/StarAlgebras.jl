# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# Example implementation of a Gram matrix that relies on `StarAlgebras.QuadraticForm` for interacting with
# `StarAlgebras.AlgebraElement`s.
struct Gram{T,B}
    matrix::T
    basis::B
end

function Base.eltype(g::Gram)
    T = eltype(g.matrix)
    basis_eltype = eltype(basis(g))
    U = if basis_eltype <: SA.AlgebraElement
        # We will multiply with the coefficients of these `AlgebraElement`
        promote_type(T, eltype(basis_eltype))
    else
        # We will multiply with the basis elements which will be keys of
        # the `SparseCoefficients` so we won't multiply with any other coefficient
        T
    end
    return MA.promote_operation(+, U, U)
end

SA.basis(g::Gram) = g.basis
Base.getindex(g::Gram, i, j) = g.matrix[i, j]
