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

@testset "Promotion of univariate to bivariate" begin
    alg_x = univariate_algebra(true)
    o = one(alg_x)
    @test isone(o)
    @test SA.coeffs(o).isless == isless
    a = SA.AlgebraElement(
        SA.SparseCoefficients(
            collect(Iterators.take(SA.object(SA.basis(alg_x)), 2)),
            [1, 5],
            isless,
        ),
        alg_x,
    )

    alg_y = univariate_algebra(false)
    b = SA.AlgebraElement(
        SA.SparseCoefficients(
            collect(Iterators.take(SA.object(SA.basis(alg_y)), 2)),
            [2, -1],
            isless,
        ),
        alg_y,
    )

    bi_alg = bivariate_algebra()
    c = SA.AlgebraElement(
        SA.SparseCoefficients(
            collect(Iterators.take(SA.object(SA.basis(bi_alg)), 3)),
            [3, -1, 5],
            grlex,
        ),
        bi_alg,
    )

    bas_x = SA.basis(alg_x)
    bas_y = SA.basis(alg_y)
    bas_bi = SA.basis(bi_alg)

    err = ErrorException("Bases $bas_x and $bas_bi are different and do not support promotion.")
    @test_throws err SA.promote_basis(bas_x, bas_bi)

    @test c == a + b

    bi_a, _b = SA.promote_basis(a, alg_y)
    @test _b == bi_alg
    @test parent(bi_a) == _b
    bi_b, _b = SA.promote_basis(b, SA.basis(alg_x))
    @test parent(bi_b) == bi_alg
    @test SA.basis(bi_b) == _b
    @test parent(bi_a) == parent(bi_b)

    @test c == bi_a + bi_b

    @testset "SubBasis" begin
        sub_x = SA.SubBasis(bas_x, [1, 3])
        sub_y = SA.SubBasis(bas_y, [2, 3])
        _sub_x, _y = SA.promote_basis(sub_x, bas_y)
        @test parent(_sub_x) == bas_bi
        @test _sub_x.keys == [(1, 0), (3, 0)]
        @test _sub_x.is_sorted
        @test _y == bas_bi
        _sub_x, _sub_y = SA.promote_basis(sub_x, sub_y)
        @test parent(_sub_x) == bas_bi
        @test _sub_x.is_sorted
        @test _sub_x.keys == [(1, 0), (3, 0)]
        @test parent(_sub_y) == bas_bi
        @test _sub_y.is_sorted
        @test _sub_y.keys == [(0, 2), (0, 3)]
    end
end
