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
            for (α, β) in ((2, 0), (0, 0), (1, 1), (2, 3), (0, 1), (0, 3))
                output = iszero(β) ? similar(result) : fill(f, size(result))
                @test (@inferred LinearAlgebra.mul!(output, A, b, α, β)) ===
                      output
                @test [[SA.coeffs(r)[i] for i in 1:5] for r in output] == [α .* c + β .* [1, 2, 0, 0, 0] for c in expected]
                @test (A, b) == originals
            end
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
        @test (@inferred LinearAlgebra.mul!(
            output,
            reshape([a, b], 2, 1),
            B,
            2,
            3,
        )) === output
        @test SA.coeffs(output[1]) == [0, 0, 40, 0, 0]
        @test SA.coeffs(output[2]) == [0, 0, 0, 60, 0]
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

@testset "Mixed sparse and dense matrix products" begin
    alg = SA.StarAlgebra(
        Monomial((0, 0)),
        SA.FixedBasis([Monomial((i, 0)) for i in 0:4]),
    )
    N = sparse([0 2; 3 0])
    for T in (Int, BigInt), storage in (identity, sparse, sparse_coefficients)
        f = SA.AlgebraElement(storage(T[1, 2, 0, 0, 0]), alg)
        g = SA.AlgebraElement(storage(T[0, 1, 1, 0, 0]), alg)
        originals = MA.mutable_copy.((f, g))
        P = SparseMatrixCSC(2, 2, Int32[1, 2, 3], Int32[2, 1], [f, g])
        identity_matrix =
            SparseMatrixCSC(2, 2, Int32[1, 2, 3], Int32[1, 2], [1, 1])
        prepared = @inferred identity_matrix * P
        @test typeof(prepared) === typeof(P)
        @test nonzeros(prepared) == [f, g]
        @test nonzeros(prepared) !== nonzeros(P)
        @test rowvals(prepared) == Int32[2, 1]
        @test rowvals(prepared) !== rowvals(P)
        @test prepared.colptr == Int32[1, 2, 3]
        @test prepared.colptr !== P.colptr
        for (A, B, expected) in (
            (N, [f, g], [2g, 3f]),
            (P, [2, 3], [3g, 2f]),
            (P, [f, g], [g * g, f * f]),
            (N, [f g; g f], [2g 2f; 3f 3g]),
            ([f g; g f], N, [3g 2f; 3f 2g]),
            (P, [2 3; 4 5], [4g 5g; 2f 3f]),
            ([2 3; 4 5], P, [3f 2g; 5f 4g]),
        )
            @test (@inferred A * B) == expected
            @test (@inferred MA.operate(*, A, B)) == expected
            output = similar(expected)
            @test (@inferred LinearAlgebra.mul!(output, A, B)) === output
            @test output == expected
            output = fill(f, size(expected))
            @test (@inferred LinearAlgebra.mul!(output, A, B, 2, 3)) === output
            @test output == [2p + 3f for p in expected]
            @test (f, g) == originals
        end
        @test +f === f
        @test *(f) === f
        t = SA.Term(alg, 2, one(T))
        @test +t === t
        @test *(t) === t

        empty = SparseMatrixCSC(2, 2, Int32[1, 1, 1], Int32[], typeof(f)[])
        @test (@inferred empty * [f, g]) == [zero(f), zero(f)]
        @test (@inferred [f g] * empty) == [zero(f) zero(f)]
        explicit = SparseMatrixCSC(2, 2, [1, 2, 2], [2], [zero(f)])
        @test (@inferred explicit * [2, 3]) == [zero(f), zero(f)]
        @test nnz(@inferred sparse([1 0; 0 1]) * explicit) == 1
        @test isempty(@inferred zeros(Int, 0, 2) * P)
    end

    a, b = SA.Term(alg, 2, 2), SA.Term(alg, 3, 3)
    P = SparseMatrixCSC(2, 2, [1, 2, 3], [2, 1], [a, b])
    output = Matrix{typeof(zero(Int, alg))}(undef, 2, 2)
    @test (@inferred LinearAlgebra.mul!(output, [1 2; 3 4], P)) === output
    @test output == [
        SA.algebra_element(2a) SA.algebra_element(b);
        SA.algebra_element(4a) SA.algebra_element(3b)
    ]
end

@testset "Sparse matrix products" begin
    alg = SA.StarAlgebra(
        Monomial((0, 0)),
        SA.FixedBasis([Monomial((i, 0)) for i in 0:4]),
    )
    N = SparseMatrixCSC(2, 2, Int32[1, 2, 3], Int32[2, 1], [3, 2])
    for T in (Int, BigInt), storage in (identity, sparse, sparse_coefficients)
        f = SA.AlgebraElement(storage(T[1, 2, 0, 0, 0]), alg)
        g = SA.AlgebraElement(storage(T[0, 1, 1, 0, 0]), alg)
        originals = MA.mutable_copy.((f, g))
        P = SparseMatrixCSC(2, 2, Int32[1, 2, 3], Int32[2, 1], [f, g])
        for (A, B, expected) in
            ((N, P, [2f, 3g]), (P, N, [3g, 2f]), (P, P, [g * f, f * g]))
            result = @inferred A * B
            @test result isa SparseMatrixCSC{typeof(f),Int32}
            @test result.colptr == Int32[1, 2, 3]
            @test rowvals(result) == Int32[1, 2]
            @test nonzeros(result) == expected
            @test nonzeros(@inferred MA.operate(*, A, B)) == expected
            output = copy(P)
            @test (@inferred LinearAlgebra.mul!(output, A, B)) === output
            @test rowvals(output) == Int32[1, 2]
            @test nonzeros(output) == expected
            @test (@inferred MA.operate_to!(output, *, A, B, 2)) === output
            @test nonzeros(output) == [2p for p in expected]
            dense = Matrix{typeof(f)}(undef, 2, 2)
            @test (@inferred LinearAlgebra.mul!(dense, A, B)) === dense
            @test dense == [expected[1] zero(f); zero(f) expected[2]]
            @test (f, g) == originals
        end
        empty = spzeros(typeof(f), Int32, 2, 2)
        @test nnz(@inferred empty * P) == 0
        @test nnz(@inferred P * empty) == 0
        @test nnz(@inferred empty * empty) == 0
        @test !SA._matrix_product_mightalias(copy(empty), transpose(empty))
        @test_throws ArgumentError LinearAlgebra.mul!(P, P, P)
        aliased =
            SparseMatrixCSC(2, 2, P.colptr, copy(rowvals(P)), copy(nonzeros(P)))
        @test_throws ArgumentError LinearAlgebra.mul!(aliased, P, P)
        aliased = SparseMatrixCSC(
            2,
            2,
            copy(empty.colptr),
            rowvals(empty),
            copy(nonzeros(empty)),
        )
        @test_throws ArgumentError LinearAlgebra.mul!(aliased, empty, P)
        @test nonzeros(P) == [f, g]
    end
end

@testset "Matrix products preserve factor order" begin
    basis = SA.DiracBasis(["", "a", "b", "ab", "ba"])
    alg = SA.StarAlgebra("", SA.DiracMStructure(basis, *))
    A, B = [1 2; 0 1], [1 0; 3 1]
    a = SA.algebra_element(SA.Term(alg, "a", A))
    b = SA.algebra_element(SA.Term(alg, "b", B))
    scalar = SA.algebra_element(SA.Term(alg, "a", 2))
    @test "b" * scalar == SA.algebra_element(SA.Term(alg, "ba", 2))
    @test scalar * "b" == SA.algebra_element(SA.Term(alg, "ab", 2))
    result = @inferred MA.operate(*, [a b], [b, a])
    @test keys(SA.coeffs(only(result))) == ["ab", "ba"]
    @test values(SA.coeffs(only(result))) == [A * B, B * A]
    left = SparseMatrixCSC(1, 2, [1, 2, 3], [1, 1], [a, b])
    @test (@inferred MA.operate(*, left, [b, a])) == result
    output = [SA.algebra_element(SA.Term(alg, "", A))]
    @test (@inferred LinearAlgebra.mul!(output, [a b], [b, a], true, true)) ===
          output
    @test keys(SA.coeffs(only(output))) == ["", "ab", "ba"]
    @test values(SA.coeffs(only(output))) == [A, A * B, B * A]
    result = @inferred MA.operate(*, [a b], reshape([b, a], 2, 1))
    @test keys(SA.coeffs(only(result))) == ["ab", "ba"]
    @test values(SA.coeffs(only(result))) == [A * B, B * A]
    right = SparseMatrixCSC(2, 1, [1, 3], [1, 2], [b, a])
    @test (@inferred MA.operate(*, [a b], right)) == result
    product = @inferred left * right
    @test only(nonzeros(product)) == only(result)
    output = Matrix{typeof(a)}(undef, 1, 1)
    @test (@inferred LinearAlgebra.mul!(output, left, right)) === output
    @test output == result
    @test A == [1 2; 0 1]
    @test B == [1 0; 3 1]
end

end
