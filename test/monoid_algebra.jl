# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

@testset "Free monoid algebra" begin
    alph = [:a, :b, :c]
    A★ = FreeWords(alph)
    B = SA.DiracBasis(A★)
    RG = StarAlgebra(A★, B)
    @test basis(RG) === B

    # no caching
    fB = SA.FixedBasis(basis(RG); n = nwords(A★, 0, 8))

    fRG = StarAlgebra(A★, fB)

    g = one(fRG)
    @test isone(g)

    @test one(fRG) == g
    @test iszero(zero(fRG))
    @test zero(g) == zero(fRG)
    @test iszero(0 * g)
    @test isone(*(g, g, g))

    @testset "Translations between bases" begin
        Z = zero(RG)
        fZ = zero(fRG)

        @test coeffs(fZ) isa SparseVector

        @test coeffs(Z, fB) == coeffs(fZ)
        @test coeffs(fZ, B) == coeffs(Z)
        @test coeffs(coeffs(fZ, B), B, fB) == coeffs(fZ)
        @test coeffs(coeffs(Z, fB), fB, B) == coeffs(Z)

        let g = basis(fRG)[rand(1:121)]
            X = RG(g)
            Y = -RG(star(g))
            for i in 1:3
                coeffs(X)[basis(fRG)[rand(1:121)]] += rand(-3:3)
                coeffs(Y)[basis(fRG)[rand(1:121)]] -= rand(3:3)
            end

            @test coeffs(X) ==
                  coeffs(coeffs(X, basis(fRG)), basis(fRG), basis(RG))
            @test coeffs(Y) ==
                  coeffs(coeffs(Y, basis(fRG)), basis(fRG), basis(RG))

            fX = AlgebraElement(coeffs(X, basis(fRG)), fRG)
            fY = AlgebraElement(coeffs(Y, basis(fRG)), fRG)

            @test dot(X, Y) == dot(fX, fY) == dot(coeffs(X), coeffs(Y))
            @test dot(fX, coeffs(fY)) == dot(coeffs(fX), fY)

            @test dot(X, X) ≈ norm(X)^2 ≈ dot(coeffs(X), coeffs(X))
            @test dot(fX, fX) ≈ norm(fX)^2 ≈ dot(coeffs(fX), coeffs(fX))

            @test coeffs(fX) ==
                  coeffs(coeffs(fX, basis(RG)), basis(RG), basis(fRG))
            @test coeffs(fY) ==
                  coeffs(coeffs(fY, basis(RG)), basis(RG), basis(fRG))

            @test coeffs(X * Y, basis(fRG)) == coeffs(fX * fY)
            @test coeffs(X * X, basis(fRG)) ==
                  coeffs(fX^2) ==
                  coeffs(X^2, basis(fRG))

            @test coeffs(Y * Y, basis(fRG)) ==
                  coeffs(fY^2) ==
                  coeffs(Y^2, basis(fRG))

            MA.operate_to!(Z, *, X, Y)
            MA.operate_to!(fZ, *, fX, fY)

            @test coeffs(Z) == coeffs(fZ, basis(RG))
            @test coeffs(fZ) == coeffs(Z, basis(fRG))

            @test coeffs(2 * X * Y) == coeffs(MA.operate_to!(Z, *, 2, Z))
            @test coeffs(2 * fX * fY) == coeffs(MA.operate_to!(fZ, *, 2, fZ))
        end
    end
    @testset "mutable arithmetic" begin
        a = let l = 12, R = 7
            support = [Word(alph, rand(1:length(alph), rand(0:R))) for _ in 1:l]
            vals = rand(-3:3, l)
            c = SA.SparseCoefficients(support, vals)
            MA.operate!(SA.canonical, c)
            AlgebraElement(c, RG)
        end
        b = let l = 7, R = 3
            support = [Word(alph, rand(1:length(alph), rand(0:R))) for _ in 1:l]
            vals = rand(-3:3, l)
            c = SA.SparseCoefficients(support, vals)
            MA.operate!(SA.canonical, c)
            AlgebraElement(c, RG)
        end

        let d = deepcopy(a)
            MA.operate!(zero, d)
            MA.operate_to!(d, -, a)

            d = deepcopy(a)
            @test !iszero(d)
            @test @allocated(MA.operate!(zero, d)) == 0
            @test iszero(d)

            @test @allocated(MA.operate_to!(d, -, a)) == 0
            @test d == -a
        end

        let d = zero(a)
            MA.operate_to!(d, +, a, b)
            @test d == a + b
            MA.operate_to!(d, +, d, b)
            @test d == a + 2b
            MA.operate_to!(d, +, b, d)
            @test d == a + 3b
            MA.operate_to!(d, +, a, b)
            @test d == a + b
        end

        let d = deepcopy(a)
            MA.operate_to!(d, *, 2, a)
            @test d == 2a
            MA.operate_to!(d, *, 2, d)
            @test d == 4a
            MA.operate_to!(d, *, a, b)
            @test d == a * b
            @test_throws ArgumentError MA.operate_to!(d, *, d, b)

            d = deepcopy(a)
            MA.operate_to!(d, *, 2, a) # preallocates memory in coefficients
            @test @allocated(MA.operate_to!(d, *, 2, a)) == 0
            @test d == 2a

            @test @allocated(MA.operate_to!(d, *, 2, d)) == 0
            @test d == 4a

            MA.operate_to!(d, *, a, b)
            @test d == a * b
            MA.operate!(SA.UnsafeAddMul(SA.mstructure(RG)), d, a, b)
            @test d == 2 * a * b

            MA.operate_to!(d, *, a, b)
            @test d == a * b
            d = a + b
            MA.operate!(SA.UnsafeAddMul(SA.mstructure(RG)), d, a, b, 3)
            @test d == a + b + 3 * a * b
        end
    end
end
