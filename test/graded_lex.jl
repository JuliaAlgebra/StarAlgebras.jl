# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# Example with Graded Lex Ordering
using Test
import StarAlgebras as SA

@testset "Graded Lex" begin
    alg = bivariate_algebra()
    o = one(alg)
    @test isone(o)
    @test SA.coeffs(o).isless == grlex
    a = SA.AlgebraElement(
        SA.SparseCoefficients(
            collect(Iterators.take(SA.object(SA.basis(alg)), 3)),
            [2, -1, 3],
            grlex,
        ),
        alg,
    )
    c = a * a
    @test c.coeffs.values == [4, -4, 12, 1, -6, 9]
    c_keys = [(0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0)]
    @test c.coeffs.basis_elements == c_keys
    sub = SA.sub_basis(SA.basis(alg), c_keys)
    @test sub.is_sorted
    for (i, key) in enumerate(c_keys)
        el = SA.basis(alg)[key]
        @test SA.key_index(sub, key) == i
        @test sub[el] == i
        @test sub[i] == el
    end
end
