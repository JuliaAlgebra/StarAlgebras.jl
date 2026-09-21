# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module TestLinearAlgebra

using Test
using SparseArrays
import LinearAlgebra
import StarAlgebras as SA
import MutableArithmetics as MA

include(joinpath(@__DIR__, "..", "examples", "bivariate.jl"))
SA.star(m::Monomial) = m

function sparse_coefficients(v)
    return SA.SparseCoefficients(findall(!iszero, v), filter(!iszero, v))
end

@testset "Matrix products retain their parent" begin
    alg = SA.StarAlgebra(
        Monomial((0, 0)),
        SA.FixedBasis([Monomial((i, 0)) for i in 0:4]),
    )
    for T in (Int, BigInt), storage in (identity, sparse, sparse_coefficients)
        f = SA.AlgebraElement(storage(T[1, 2, 0, 0, 0]), alg)
        g = SA.AlgebraElement(storage(T[0, 1, 1, 0, 0]), alg)
        for (A, b, expected) in (
            ([2 -1; 1 3], [f, g], [[2, 3, -1, 0, 0], [1, 5, 3, 0, 0]]),
            ([f g; g f], [2, 3], [[2, 7, 3, 0, 0], [3, 8, 2, 0, 0]]),
            ([f g; g f], [f, g], [[1, 4, 5, 2, 1], [0, 2, 6, 4, 0]]),
            (
                [2 -1; 1 3],
                [f g; g f],
                reshape(
                    [
                        [2, 3, -1, 0, 0],
                        [1, 5, 3, 0, 0],
                        [-1, 0, 2, 0, 0],
                        [3, 7, 1, 0, 0],
                    ],
                    2,
                    2,
                ),
            ),
            (
                [f g; g f],
                [2 3; 3 1],
                reshape(
                    [
                        [2, 7, 3, 0, 0],
                        [3, 8, 2, 0, 0],
                        [3, 7, 1, 0, 0],
                        [1, 5, 3, 0, 0],
                    ],
                    2,
                    2,
                ),
            ),
            (
                [f g; g f],
                [f g; g f],
                reshape(
                    [
                        [1, 4, 5, 2, 1],
                        [0, 2, 6, 4, 0],
                        [0, 2, 6, 4, 0],
                        [1, 4, 5, 2, 1],
                    ],
                    2,
                    2,
                ),
            ),
        )
            originals = MA.mutable_copy.((A, b))
            result = @inferred MA.operate(*, A, b)
            @test [[SA.coeffs(r)[i] for i in 1:5] for r in result] == expected
            @test all(parent(r) === alg for r in result)
            @test (@inferred A * b) == result
            @test (A, b) == originals
            output = similar(result)
            @test (@inferred LinearAlgebra.mul!(output, A, b)) === output
            @test output == result
            second = MA.mutable_copy(result[2])
            MA.operate!(zero, result[1])
            @test result[2] == second
            @test (A, b) == originals
        end
        @test isempty(@inferred MA.operate(*, zeros(Int, 0, 2), [f, g]))
    end

    f = SA.AlgebraElement(SA.SparseCoefficients([1, 2], [1, 2]), alg)
    g = SA.AlgebraElement([0, 1, 1, 0, 0], alg)
    result = @inferred MA.operate(*, [f f; f f], [g, g])
    @test SA.coeffs(result[1]) == [0, 2, 6, 4, 0]
    @test result[1] == result[2]

    a, b = SA.Term(alg, 2, 2), SA.Term(alg, 3, 3)
    for shape in ((2,), (2, 1))
        output = Array{typeof(zero(Int, alg))}(undef, shape)
        B = fill(SA.Term(alg, 2, 4), ntuple(_ -> 1, length(shape)))
        @test (@inferred MA.operate_to!(
            output,
            *,
            reshape([a, b], 2, 1),
            B,
        )) === output
        @test SA.coeffs(output[1]) == [0, 0, 8, 0, 0]
        @test SA.coeffs(output[2]) == [0, 0, 0, 12, 0]
    end

    inputs = [f, copy(f)]
    originals = MA.mutable_copy(inputs)
    @test_throws ArgumentError LinearAlgebra.mul!(inputs, [1 2; 3 4], inputs)
    @test inputs == originals
    @test_throws ArgumentError LinearAlgebra.mul!(
        inputs,
        reshape(inputs, 2, 1),
        [1],
    )
    @test inputs == originals
    output = [copy(f)]
    @test_throws DimensionMismatch LinearAlgebra.mul!(
        output,
        [1 2; 3 4],
        inputs,
    )
    @test output == [f]
    @test_throws DimensionMismatch LinearAlgebra.mul!(
        output,
        [1 2],
        inputs[1:1],
    )
    @test output == [f]
end

@testset "Matrix products preserve factor order" begin
    basis = SA.DiracBasis(["", "a", "b", "ab", "ba"])
    alg = SA.StarAlgebra("", SA.DiracMStructure(basis, *))
    A, B = [1 2; 0 1], [1 0; 3 1]
    a = SA.algebra_element(SA.Term(alg, "a", A))
    b = SA.algebra_element(SA.Term(alg, "b", B))
    result = @inferred MA.operate(*, [a b], [b, a])
    @test keys(SA.coeffs(only(result))) == ["ab", "ba"]
    @test values(SA.coeffs(only(result))) == [A * B, B * A]
    result = @inferred MA.operate(*, [a b], reshape([b, a], 2, 1))
    @test keys(SA.coeffs(only(result))) == ["ab", "ba"]
    @test values(SA.coeffs(only(result))) == [A * B, B * A]
    @test A == [1 2; 0 1]
    @test B == [1 0; 3 1]
end

end
