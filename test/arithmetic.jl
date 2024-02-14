@testset "Arithmetic" begin
    G = PermGroup(perm"(1,2,3)", perm"(1,2)")
    b = StarAlgebras.Basis{UInt8}(collect(G))
    l = length(b)
    RG = StarAlgebra(G, b, (l, l))

    @test_broken contains(sprint(show, RG), "*-algebra of Permutation group")

    @testset "Module structure" begin
        a = AlgebraElement(ones(Int, order(G)), RG)

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
        a = AlgebraElement(ones(Int, order(G)), RG)
        b = sum(sign(g) * RG(g) for g in G)

        @test a ==
              AlgebraElement(ones(Int, order(G)), RG) ==
              sum(RG(g) for g in G)

        @test 1 / 2 * (coeffs(a + b)) == [1.0, 0.0, 1.0, 0.0, 1.0, 0.0]

        g, h = gens(G)
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
        @test LinearAlgebra.dot(a, coeffs(a)) ≈
              norm(a)^2 ≈
              LinearAlgebra.dot(coeffs(a), a)

        @test a * b == StarAlgebras.mul!(a, a, b)

        @test aug(a) == 3
        @test aug(b) == -1
        @test aug(a) * aug(b) == aug(a * b) == aug(b * a)

        z = sum((one(RG) - RG(g)) * star(one(RG) - RG(g)) for g in G)
        @test aug(z) == 0

        @test supp(z) == StarAlgebras.basis(parent(z))
        @test supp(RG(1) + RG(g)) == [one(G), g]
        @test supp(a) == [one(G), g, h]

        @testset "Projections in star algebras" begin
            b = StarAlgebras.basis(RG)
            l = length(b)
            P = sum(RG(g) for g in b) // l
            @test P * P == P

            P3 = 2 * sum(RG(g) for g in b if sign(g) > 0) // l
            @test P3 * P3 == P3

            PAlt = sum(sign(g) * RG(g) for g in b) // l
            @test PAlt * PAlt == PAlt

            @test P3 * PAlt == PAlt * P3

            P2 = (RG(1) + RG(b[2])) // 2
            @test P2 * P2 == P2

            @test P2 * P3 == P3 * P2 == P

            P2m = (RG(1) - RG(b[2])) // 2
            @test P2m * P2m == P2m

            @test P2m * P3 == P3 * P2m == PAlt
            @test iszero(P2m * P2)
        end
    end

    @testset "Mutable API and trivial mstructure" begin
        A = [:a, :b, :c]
        b = StarAlgebras.Basis{UInt16}(words(A, radius=8))
        l = findfirst(w -> length(w) > 4, b) - 1

        RG = StarAlgebra(one(first(b)), b)

        @test basis(RG) === b
        @test basis(RG.mstructure) === basis(RG)

        RGc = StarAlgebra(one(first(b)), b, (l, l))
        @test basis(RGc) === b
        @test basis(RGc.mstructure) === basis(RGc)

        @test all(RG.mstructure[1:121, 1:121] .== RGc.mstructure)

        Z = zero(RG)
        W = zero(RGc)

        let g = b[rand(1:121)]
            X = RG(g)
            Y = -RG(star(g))
            for i in 1:3
                X[b[rand(1:121)]] += rand(-3:3)
                Y[b[rand(1:121)]] -= rand(3:3)
            end

            Xc = AlgebraElement(coeffs(X), RGc)
            Yc = AlgebraElement(coeffs(Y), RGc)

            @test coeffs(X * Y) ==
                  coeffs(Xc * Yc) ==
                  coeffs(StarAlgebras.mul!(Z, X, Y))

            @test coeffs(X^2) == coeffs(Xc^2) == coeffs(X * X)
            @test coeffs(Y^2) == coeffs(Yc^2) == coeffs(Y * Y)

            @test coeffs(Z) == StarAlgebras.mul!(
                coeffs(W),
                coeffs(X),
                coeffs(Y),
                RG.mstructure,
            )
            @test coeffs(Z) == coeffs(W)
            @test coeffs(Z) == StarAlgebras.mul!(
                coeffs(W),
                coeffs(X),
                coeffs(Y),
                RGc.mstructure,
            )
            @test coeffs(Z) == coeffs(W)

            StarAlgebras.zero!(W)
            StarAlgebras.fmac!(coeffs(W), coeffs(X), coeffs(Y), RG.mstructure)

            @test coeffs(2 * X * Y) == coeffs(StarAlgebras.mul!(W, W, 2))
        end
    end

    @testset "mutable arithmetic" begin
        A = [:a, :b, :c]
        bas = StarAlgebras.Basis{UInt16}(words(A, radius=4))
        l = findfirst(w -> length(w) > 2, bas) - 1
        RG = StarAlgebra(one(first(bas)), bas, (l, l))

        a = let c = rand(-3:3, l)
            resize!(c, length(bas))
            c[l:end] .= 0
            AlgebraElement(c, RG)
        end
        b = let c = rand(-3:3, l)
            resize!(c, length(bas))
            c[l:end] .= 0
            AlgebraElement(c, RG)
        end

        let d = deepcopy(a)
            StarAlgebras.zero!(d)
            StarAlgebras.neg!(d, a)

            d = deepcopy(a)
            @test !iszero(d)
            @test @allocated(StarAlgebras.zero!(d)) == 0
            @test iszero(d)

            @test @allocated(StarAlgebras.neg!(d, a)) == 0
            @test d == -a
        end

        let d = deepcopy(a)
            StarAlgebras.add!(d, d, b)
            StarAlgebras.add!(d, b, d)
            StarAlgebras.add!(d, a, b)

            d = deepcopy(a)
            @test @allocated(StarAlgebras.add!(d, d, b)) == 0
            @test d == a + b

            d = deepcopy(a)
            @test @allocated(StarAlgebras.add!(d, b, d)) == 0
            @test d == a + b

            @test @allocated(StarAlgebras.add!(d, a, b)) == 0
            @test d == a + b
        end

        let d = deepcopy(a)
            StarAlgebras.mul!(d, d, 2)
            StarAlgebras.mul!(d, a, 2)
            StarAlgebras.mul!(d, a, b)
            d = deepcopy(a)
            StarAlgebras.mul!(d, d, b)

            d = deepcopy(a)
            @test @allocated(StarAlgebras.mul!(d, d, 2)) == 0
            @test d == 2a

            @test @allocated(StarAlgebras.mul!(d, a, 2)) == 0
            @test d == 2a

            @test @allocated(StarAlgebras.mul!(d, a, b)) == 32
            @test d == a * b

            d = deepcopy(a)
            @test @allocated(StarAlgebras.mul!(d, d, b)) != 0
            z = StarAlgebras.mul!(d, d, b)
            @test z == a * b
            @test z !== d
        end
    end
end

@testset "Group Algebra caching" begin
    A = [:a, :b, :c]
    b = StarAlgebras.Basis{UInt8}(words(A, radius=4))
    k = findfirst(w -> length(w) == 3, b) - 1

    RG = StarAlgebra(Word(A, Int[]), b, (k, k))
    @test RG isa StarAlgebra

    D = sum(RG(b[i]) for i in 1:k)
    @test D isa AlgebraElement
    g = one(RG)
    @test isone(g)

    @test one(RG) == g
    @test iszero(zero(RG))
    @test 0 * g == zero(RG)
    @test iszero(0 * g)

    h = RG(b[3])

    @test D * one(RG) == D
    @test one(RG) * D == D

    @test supp(D) == b[1:k]

    @test_throws StarAlgebras.ProductNotDefined all(!iszero, RG.mstructure.table)

    @test D * D isa AlgebraElement

    @test all(!iszero, RG.mstructure.table)

    RG = StarAlgebra(Word(A, Int[]), b, (k, k), precompute=true)
    @test all(!iszero, RG.mstructure.table)
end
