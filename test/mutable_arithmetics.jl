# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module TestMutableArithmetics

using Test
using SparseArrays
using StarAlgebras
import StarAlgebras as SA
import MutableArithmetics as MA

include(joinpath(@__DIR__, "..", "examples", "acoeffs.jl"))

@testset "Mutable copies and operation routing" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0, 4.0]))
    for c in (
        BigInt[1, 0, 3, 0],
        sparsevec([1, 3], BigInt[1, 3], 4),
        SA.SparseCoefficients([1, 3], BigInt[1, 3]),
    )
        a = AlgebraElement(c, alg)
        @test (@inferred MA.promote_operation(zero, typeof(a))) ===
              typeof(zero(a))
        original = deepcopy(a)
        b = MA.copy_if_mutable(a)
        @test b !== a
        @test parent(b) === parent(a)
        @test SA.coeffs(b) !== c
        @test SA.coeffs(b)[1] !== c[1]
        MA.operate!(+, SA.coeffs(b)[1], 1)
        @test SA.coeffs(b)[1] == 2
        @test c[1] == 1

        @test MA.operate!!(remove_leading_term, b) === b
        @test iszero(SA.coeffs(b)[3])
        @test SA.coeffs(b)[1] == 2
        @test MA.operate!!(zero, b) === b
        @test iszero(b)
        @test a == original

        # Missing in-place implementations must surface instead of allocating.
        @test MA.mutability(typeof(a), +, typeof(a), typeof(a)) isa
              MA.IsMutable
        @test_throws ErrorException MA.operate!!(+, a, a)
        @test_throws ErrorException MA.operate!!(MA.sub_mul, a, 2, a)
        @test a == original
    end
end

@testset "Zero promotion with tuple storage" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0, 4.0]))
    a = AlgebraElement(SA.SparseCoefficients((1, 3), (big(1), big(3))), alg)
    z = zero(a)
    @test (@inferred MA.promote_operation(zero, typeof(a))) === typeof(z)
    @test typeof(z) !== typeof(a)
    @test parent(z) === alg

    r = @inferred MA.operate!!(zero, a)
    @test iszero(r)
    @test typeof(r) === typeof(z)
    @test parent(r) === alg
    @test keys(SA.coeffs(a)) == (1, 3)
    @test values(SA.coeffs(a)) == (1, 3)
end

@testset "Coefficient storage mutation routing" begin
    v = sparsevec([1], BigInt[2], 3)
    @test MA.copy_if_mutable(v) === v

    c = SA.SparseCoefficients([1, 3], BigInt[2, 3])
    b = MA.copy_if_mutable(c)
    @test b !== c
    MA.operate!(+, b[1], 1)
    @test c[1] == 2
    result = MA.operate!!(+, b, c)
    @test result === b
    @test result[1] == 5
    @test c[1] == 2
    @test MA.operate!!(remove_leading_term, b) === b
    @test iszero(b[3])
    @test MA.operate!!(zero, b) === b
    @test iszero(b)
end

@testset "Sparse coefficient merge" begin
    for T in (Int, BigInt), alias in (:none, :left, :right)
        x = SA.SparseCoefficients([1, 3, 5], T[2, 3, 4])
        y = SA.SparseCoefficients([2, 3, 6], T[5, -3, 7])
        out = alias === :left ? x : alias === :right ? y : zero(x)
        @test MA.operate_to!(out, +, x, y) === out
        @test keys(out) == [1, 2, 5, 6]
        @test values(out) == [2, 5, 4, 7]
        alias === :left || @test values(x) == [2, 3, 4]
        alias === :right || @test values(y) == [5, -3, 7]
        @test MA.operate!!(+, out, out) === out
        @test values(out) == [4, 10, 8, 14]
        @test MA.operate!!(+, out, -out) === out
        @test isempty(keys(out))
        @test isempty(values(out))
    end

    x = SA.SparseCoefficients([3, 1, 3], [1, 2, -1])
    y = SA.SparseCoefficients([2, 1], [3, -2], >)
    @test MA.operate!!(+, x, y) === x
    @test keys(x) == [2]
    @test values(x) == [3]
    @test keys(y) == [2, 1]
    @test values(y) == [3, -2]

    x = SA.SparseCoefficients([1], [1])
    y = SA.SparseCoefficients([1, 2], [0.5, 2.5])
    out = MA.operate!!(+, x, y)
    @test out !== x
    @test values(out) == [1.5, 2.5]
    @test values(x) == [1]
end

@testset "Custom coefficient storage" begin
    alg = StarAlgebra(1.0, SA.FixedBasis{Float64,UInt32}([1.0, 2.0, 3.0, 4.0]))
    a = AlgebraElement(ACoeffs(BigInt[1, 0, 3, 0]), alg)
    b = MA.copy_if_mutable(a)
    @test typeof(b) === typeof(a)
    @test parent(b) === parent(a)
    @test SA.coeffs(b).vals !== SA.coeffs(a).vals
    MA.operate!(+, SA.coeffs(b)[1], 1)
    @test SA.coeffs(b)[1] == 2
    @test SA.coeffs(a)[1] == 1
    @test MA.operate!!(remove_leading_term, b) === b
    @test SA.coeffs(b).vals == [2, 0, 0, 0]
    @test MA.operate!!(zero, b) === b
    @test iszero(b)
    @test SA.coeffs(a).vals == [1, 0, 3, 0]

    # Allocating negation works, but missing MA support must not be hidden.
    @test (-SA.coeffs(a)).vals == [-1, 0, -3, 0]
    @test_throws MethodError MA.operate!!(-, SA.coeffs(a))
end

@testset "Coefficient mapping" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0, 4.0]))
    for T in (Int, BigInt), alias in (false, true)
        for c in (
            T[2, 0, 3, 0],
            sparsevec([1, 3], T[2, 3], 4),
            SA.SparseCoefficients([1, 3], T[2, 3]),
        )
            input = AlgebraElement(c, alg)
            output = alias ? input : copy(input)
            @test @inferred(
                SA.map_coefficients_to!(output, c -> c - 2, input)
            ) === output
            @test [SA.coeffs(output)[i] for i in 1:4] == [0, 0, 1, 0]
            if !alias
                @test [SA.coeffs(input)[i] for i in 1:4] == [2, 0, 3, 0]
            end
            @test SA.map_coefficients!(c -> 2c, output; nonzero = true) ===
                  output
            @test [SA.coeffs(output)[i] for i in 1:4] == [0, 0, 2, 0]
            @test SA.map_coefficients!(zero, output) === output
            @test iszero(output)
            if c isa Union{SA.SparseCoefficients,SparseVector}
                @test isempty(collect(SA.nonzero_pairs(SA.coeffs(output))))
            end
        end
    end

    # Preserve canonical output when storage uses a different ordering.
    c = SA.SparseCoefficients([1, 2, 3], [1, 2, 3])
    out = SA.SparseCoefficients(Int[], Float64[], >)
    @test SA.map_coefficients_to!(out, c -> c / 2, c) === out
    @test keys(out) == [3, 2, 1]
    @test values(out) == [1.5, 1.0, 0.5]
    @test values(c) == [1, 2, 3]
    repeated = SA.SparseCoefficients([3, 1, 3], [1, 2, 3])
    @test SA.map_coefficients_to!(out, c -> c^2, repeated) === out
    @test keys(out) == [3, 1]
    @test values(out) == [10, 4]

    input = AlgebraElement(c, alg)
    other = AlgebraElement(copy(c), StarAlgebra(2.0, basis(alg)))
    @test_throws ArgumentError SA.map_coefficients_to!(other, zero, input)
    @test values(SA.coeffs(other)) == [1, 2, 3]
end

end
