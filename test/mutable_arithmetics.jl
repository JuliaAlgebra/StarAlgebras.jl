# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

@testset "Mutable copies and operation routing" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0, 4.0]))
    for c in (
        BigInt[1, 0, 3, 0],
        sparsevec([1, 3], BigInt[1, 3], 4),
        SA.SparseCoefficients([1, 3], BigInt[1, 3]),
    )
        a = AlgebraElement(c, alg)
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

        # Operations without an in-place implementation still allocate.
        @test MA.mutability(typeof(a), +, typeof(a), typeof(a)) isa
              MA.IsNotMutable
        @test MA.operate!!(+, a, a) == 2a
        @test MA.operate!!(MA.sub_mul, a, 2, a) == -a
        @test a == original
    end
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
    @test result[1] == 5
    @test b[1] == 3
    @test MA.operate!!(remove_leading_term, b) === b
    @test iszero(b[3])
    @test MA.operate!!(zero, b) === b
    @test iszero(b)
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
