# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module TestSum

using Test
using SparseArrays
import StarAlgebras as SA
import MutableArithmetics as MA

@testset "Sums retain their parent" begin
    alg = SA.StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0]))
    for T in (Int, BigInt)
        for c in (
            T[1, 0, 3],
            sparsevec([1, 3], T[1, 3], 3),
            SA.SparseCoefficients([1, 3], T[1, 3]),
        )
            f = SA.AlgebraElement(c, alg)
            g = 2f
            a = [f, g]
            originals = MA.mutable_copy(a)
            result = @inferred sum(a)
            @test result == 3f
            @test parent(result) === alg
            @test (@inferred MA.operate(sum, a)) == result
            @test sum(a; init = f) == 4f
            @test sum([f]) == f
            @test SA.coeffs(sum([f])) !== c
            @test a == originals

            @test sum(typeof(f)[]; init = f) == f
            @test SA.coeffs(sum(typeof(f)[]; init = f)) !== c
            matrix = [f g; g f]
            @test sum(matrix; dims = 1) == [3f 3f]
            @test sum(matrix; dims = 2) == reshape([3f, 3f], 2, 1)
            @test sum(matrix; dims = (1, 2), init = f) == fill(7f, 1, 1)
            @test a == originals
        end
    end

    terms = [SA.Term(alg, 1, 2), SA.Term(alg, 3, 4)]
    result = @inferred sum(terms)
    @test SA.coeffs(result) == [2, 0, 4]
    @test parent(result) === alg
    @test sum(terms; init = terms[1]) == SA.algebra_element(terms[1]) + result
    @test SA.coefficient.(terms) == [2, 4]

    c = SA.SparseCoefficients((1, 3), (2, 4))
    f = SA.AlgebraElement(c, alg)
    result = @inferred sum([f, f])
    @test keys(SA.coeffs(result)) == [1, 3]
    @test values(SA.coeffs(result)) == [4, 8]
    @test values(c) == (2, 4)
end

end
