@testset "Permutation group algebra" begin
    G = PermGroup(perm"(1,2,3)", perm"(1,2)")
    b = SA.DiracBasis{UInt8}(G)
    RG = StarAlgebra(G, b)

    @test contains(sprint(show, RG), "*-algebra of")

    @testset "Module structure" begin
        sc = SA.SparseCoefficients(
            collect(SA.object(RG)),
            ones(Int, length(basis(RG))),
        )
        a = AlgebraElement(sc, RG)

        @test -a isa AlgebraElement
        @test coeffs(-a) == -coeffs(a)

        @test 2 * a isa AlgebraElement
        @test eltype(2 * a) == typeof(2)
        @test coeffs(2 * a) == 2coeffs(a)

        @test eltype(2a) == Int
        y = div(2a, 2)
        @test y == a
        @test eltype(y) == Int

        @test 2.0 * a isa AlgebraElement
        @test eltype(2.0 * a) == typeof(2.0)
        @test coeffs(2.0 * a) == 2.0 * coeffs(a)

        @test coeffs(a / 2) == coeffs(a) / 2
        b = a / 2
        @test b isa AlgebraElement
        @test eltype(b) == typeof(1 / 2)
        @test coeffs(b / 2) == 0.25 * coeffs(a)

        c = a // 1

        @test eltype(c) == Rational{Int}
        @test c // 4 isa AlgebraElement
        @test c // big(4) isa AlgebraElement
        @test eltype(c // (big(4) // 1)) == Rational{BigInt}

        @test (1.0a) * 1 // 2 == (0.5a) == c // 2
    end

    @testset "Additive structure" begin
        sc = SA.SparseCoefficients(
            collect(SA.object(RG)),
            ones(Int, length(basis(RG))),
        )
        a = AlgebraElement(sc, RG)
        b = sum(sign(g) * RG(g) for g in G)

        @test a == sum(RG(g) for g in G)

        @test 1 / 2 * (coeffs(a + b)) == SA.SparseCoefficients(
            collect(SA.object(RG)),
            [1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
        )

        g, h = PermutationGroups.gens(G)
        k = g * h

        a = RG(1) + RG(g) + RG(h)
        b = RG(1) - RG(k) - RG(h)

        @test a - b == RG(g) + RG(k) + 2RG(h)

        @test a + b - 2a == b - a

        @test 1 // 2 * 2a == a
        @test a + 2a == (3 // 1) * a
        @test 2a - (1 // 1) * a == a
    end

    @testset "Multiplicative structure" begin
        for g in G, h in G
            a = RG(g)
            b = RG(h)
            @test a * b == RG(g * h)
            @test (a + b) * (a + b) == a * a + a * b + b * a + b * b
        end

        for g in G
            @test star(RG(g)) == RG(inv(g))
            @test (one(RG) - RG(g)) * star(one(RG) - RG(g)) ==
                  2 * one(RG) - RG(g) - RG(inv(g))
            @test aug(one(RG) - RG(g)) == 0
        end

        g, h = PermutationGroups.gens(G)
        k = g * h

        a = RG(1) + RG(g) + RG(h)
        b = RG(1) - RG(k) - RG(h)

        @test norm(a) == norm(b)

        @test a * b == MA.operate_to!(similar(a), *, a, b)

        @test aug(a) == 3
        @test aug(b) == -1
        @test aug(a) * aug(b) == aug(a * b) == aug(b * a)

        z = sum((one(RG) - RG(g)) * star(one(RG) - RG(g)) for g in G)
        @test aug(z) == 0

        @test supp(z) == sort(collect(basis(parent(z))))
        @test supp(RG(1) + RG(g)) == [one(G), g]
        @test supp(a) == [one(G), h, g]

        @testset "Projections in star algebras" begin
            b = basis(RG)
            l = length(b)
            P = sum(RG(g) for g in b) // l
            @test P * P == P

            P3 = 2 * sum(RG(g) for g in b if sign(g) > 0) // l
            @test P3 * P3 == P3

            PAlt = sum(sign(g) * RG(g) for g in b) // l
            @test PAlt * PAlt == PAlt

            @test P3 * PAlt == PAlt * P3

            h = Permutation(perm"(2,3)", G)
            P2 = (RG(1) + RG(h)) // 2
            @test P2 * P2 == P2

            @test P2 * P3 == P3 * P2 == P

            P2m = (RG(1) - RG(h)) // 2
            @test P2m * P2m == P2m

            @test P2m * P3 == P3 * P2m == PAlt
            @test iszero(P2m * P2)
        end
    end
end

@testset "Free monoid algebra" begin
    alph = [:a, :b, :c]
    A★ = FreeWords(alph)
    B = SA.DiracBasis{UInt16}(A★)
    RG = StarAlgebra(A★, B)
    @test basis(RG) === B

    # no caching
    fB = SA.FixedBasis(basis(RG); n = nwords(A★, 0, 8), mt = 0)
    @test fB.table.elts === fB.elts

    fRG = StarAlgebra(A★, fB)

    g = one(fRG)
    @test isone(g)

    @test one(fRG) == g
    @test iszero(zero(fRG))
    @test zero(g) == zero(fRG)
    @test_broken iszero(0 * g)

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

            @test LinearAlgebra.dot(fX, coeffs(fX)) ≈
                  norm(fX)^2 ≈
                  LinearAlgebra.dot(coeffs(fX), fX)

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

            @test coeffs(2 * X * Y) == coeffs(MA.operate_to!(Z, *, Z, 2))
            @test_broken coeffs(2 * fX * fY) ==
                         coeffs(MA.operate_to!(fZ, *, fZ, 2))
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
            MA.operate_to!(d, *, a, 2)
            @test d == 2a
            MA.operate_to!(d, *, d, 2)
            @test d == 4a
            MA.operate_to!(d, *, a, b)
            @test d == a * b
            @test_throws ArgumentError MA.operate_to!(d, *, d, b)

            d = deepcopy(a)
            @test @allocated(MA.operate_to!(d, *, d, 2)) == 0
            @test d == 2a

            @test @allocated(MA.operate_to!(d, *, a, 2)) == 0
            @test d == 2a
        end
    end

    @testset "FixedBasis caching && allocations" begin
        fB = SA.FixedBasis(B; n = nwords(A★, 8), mt = UInt32(nwords(A★, 4)))
        fRG = StarAlgebra(A★, fB)

        k = size(SA.mstructure(basis(fRG)), 1)

        y = spzeros(length(basis(fRG)))
        y[1:k] .= 1
        Y = AlgebraElement(y, fRG)
        @test_broken Y = sum(fRG(basis(fRG)[i]) for i in 1:k)

        @test Y isa AlgebraElement

        @static if VERSION ≥ v"1.9"
            star(Y)
            star(Y)
            @test (@allocations star(Y)) ≤ 4
        end

        @test supp(Y) == basis(fRG)[1:k]

        @test Y * one(fRG) == Y
        @test one(fRG) * Y == Y

        @test_throws SA.UndefRefError all(!iszero, SA.mstructure(fRG).table)

        @static if VERSION ≥ v"1.9"
            @test (@allocations Y * Y) > k^2 - 2 * k
            @test Y * Y isa AlgebraElement
            @test (@allocations Y * Y) ≤ 26
        else
            k1 = @allocated Y * Y
            @test Y * Y isa AlgebraElement
            Y * Y
            k2 = @allocated Y * Y
            @test k2 / k1 < 0.5
        end

        @test all(!iszero, SA.mstructure(fRG).table)

        @static if VERSION ≥ v"1.9"
            YY = deepcopy(Y)
            # MA.operate_to!(YY, +, Y, YY)
            # YY = deepcopy(Y)
            @test_broken @allocations(MA.operate_to!(YY, +, Y, YY)) == 0
            @test_broken YY == Y + Y

            # MA.operate_to!(YY, +, YY, Y)
            # YY = deepcopy(Y)
            @test_broken @allocations(MA.operate_to!(YY, +, YY, Y)) == 0
            @test_broken YY == Y + Y

            # MA.operate_to!(YY, +, Y, Y)
            @test_broken @allocations(MA.operate_to!(YY, +, Y, Y)) == 0
            @test_broken YY == Y + Y

            MA.operate_to!(YY, *, Y, Y)
            @test @allocations(MA.operate_to!(YY, *, Y, Y)) ≤ 25
            @test YY == Y * Y
        end
    end
end
