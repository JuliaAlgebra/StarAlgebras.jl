@testset "Permutation group algebra: arithmetic " begin
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
        @test zero(coeffs(a)) == coeffs(zero(a))

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
        @test (2a) // 2 == a
        @test a + 2a == 3.0a
        @test 2a - a // 1 == a
        @test div(3a, 2) == a

        @test coeffs(a + b - 2a) == coeffs(a) + coeffs(b) - 2coeffs(a)
        @test coeffs(2a // 2) == 2coeffs(a) // 2 == coeffs(a)
        @test coeffs(3a) == 3.0coeffs(a)
        @test coeffs(2a) - coeffs(a // 1) ==
              2coeffs(a) - coeffs(a) // 1 ==
              coeffs(a)
        @test div(3coeffs(a), 2) == coeffs(a)
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
            @test SA.aug(one(RG) - RG(g)) == 0
        end

        g, h = PermutationGroups.gens(G)
        k = g * h

        a = RG(1) + RG(g) + RG(h)
        b = RG(1) - RG(k) - RG(h)

        @test norm(a) == norm(b)

        @test a * b == MA.operate_to!(similar(a), *, a, b)

        @test SA.aug(a) == 3
        @test SA.aug(b) == -1
        @test SA.aug(a) * SA.aug(b) == SA.aug(a * b) == SA.aug(b * a)

        z = sum((one(RG) - RG(g)) * star(one(RG) - RG(g)) for g in G)
        @test SA.aug(z) == 0

        @test SA.supp(z) == sort(collect(basis(parent(z))))
        @test SA.supp(RG(1) + RG(g)) == [one(G), g]
        @test SA.supp(a) == [one(G), h, g]

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
