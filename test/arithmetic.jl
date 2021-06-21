using GroupsCore
import Random

include(joinpath(pathof(GroupsCore), "..", "..", "test", "cyclic.jl"))
StarAlgebras.star(g::GroupElement) = inv(g)

@testset "Arithmetic" begin
    G = CyclicGroup(6)
    b = StarAlgebras.Basis{UInt8}(collect(G))
    l = length(b)
    RG = StarAlgebra(G, b, (l,l))

    @testset "Module structure" begin
        a = AlgebraElement(ones(Int, order(G)), RG)

        @test -a isa AlgebraElement
        @test coeffs(-a) == -coeffs(a)

        @test 2*a isa AlgebraElement
        @test eltype(2*a) == typeof(2)
        @test coeffs(2*a) == 2coeffs(a)

        @test 2.0*a isa AlgebraElement
        @test eltype(2.0*a) == typeof(2.0)
        @test coeffs(2.0*a) == 2.0*coeffs(a)

        @test coeffs(a/2) == coeffs(a)/2
        b = a/2
        @test b isa AlgebraElement
        @test eltype(b) == typeof(1/2)
        @test coeffs(b/2) == 0.25*coeffs(a)

        c = a//1

        @test eltype(c) == Rational{Int}
        @test c//4 isa AlgebraElement
        @test c//big(4) isa AlgebraElement
        @test eltype(c//(big(4)//1)) == Rational{BigInt}

        @test (1.0a)*1//2 == (0.5a) == c//2
    end

    @testset "Additive structure" begin
        a = AlgebraElement(ones(Int, order(G)), RG)
        b = sum((-1)^isodd(g.residual)*RG(g) for g in G)

        @test a == AlgebraElement(ones(Int, order(G)), RG) == sum(RG(g) for g in G)

        @test 1/2*(a+b).coeffs == [1.0, 0.0, 1.0, 0.0, 1.0, 0.0]

        g = CyclicGroupElement(2, G)
        h = CyclicGroupElement(3, G)
        k = CyclicGroupElement(5, G)

        a = RG(1) + RG(g) + RG(h)
        b = RG(1) - RG(k) - RG(h)

        @test a - b == RG(g) + RG(k) + 2RG(h)

        @test a + b - 2a == b - a

        @test 1//2*2a == a
        @test a + 2a == (3//1)*a
        @test 2a - (1//1)*a == a
    end

    @testset "Multiplicative structure" begin

        for g in G, h in G
            a = RG(g)
            b = RG(h)
            @test a*b == RG(g*h)
            @test (a+b)*(a+b) == a*a + a*b + b*a + b*b
        end

        for g in G
            @test star(RG(g)) == RG(inv(g))
            @test (one(RG)-RG(g))*star(one(RG)-RG(g)) ==
            2*one(RG) - RG(g) - RG(inv(g))
            @test aug(one(RG) - RG(g)) == 0
        end

        g = CyclicGroupElement(2, G)
        h = CyclicGroupElement(3, G)
        k = CyclicGroupElement(5, G)

        a = RG(1) + RG(g) + RG(h)
        b = RG(1) - RG(k) - RG(h)

        # a = RG(1) + RG(perm"(2,3)") + RG(perm"(1,2,3)")
        # b = RG(1) - RG(perm"(1,2)(3)") - RG(perm"(1,2,3)")

        @test a*b == StarAlgebras.mul!(a,a,b)

        @test aug(a) == 3
        @test aug(b) == -1
        @test aug(a)*aug(b) == aug(a*b) == aug(b*a)

        z = sum((one(RG)-RG(g))*star(one(RG)-RG(g)) for g in G)
        @test aug(z) == 0

        @test supp(z) == basis(parent(z))
        @test supp(RG(1) + RG(g)) == [one(G), g]
        @test supp(a) == [one(G), g, h]

        if false
            @testset "Projections in Symm(3)" begin
                G = SymmetricGroup(3)
                b = StarAlgebras.Basis{UInt8}(collect(G))
                l = length(b)

                RG = StarAlgebra(G, b)
                @test RG isa StarAlgebra

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
    end

    @testset "Mutable API and trivial mstructure" begin

        A = [:a, :b, :c]
        b = StarAlgebras.Basis{UInt16}(words(A, radius=8))
        l = findfirst(w->length(w)>4, b)-1

        RG = StarAlgebra(one(first(b)), b)

        @test basis(RG) === b
        @test basis(RG.mstructure) === basis(RG)

        RGc = StarAlgebra(one(first(b)), b, (l, l))
        @test basis(RGc) === b
        @test basis(RGc.mstructure) === basis(RGc)

        @test all(RG.mstructure[1:121, 1:121] .== RGc.mstructure)

        Z = zero(RGc)
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

            @test coeffs(X*Y) == coeffs(Xc*Yc) == coeffs(StarAlgebras.mul!(Z, X, Y))

            @test coeffs(X^2) == coeffs(Xc^2) == coeffs(X*X)
            @test coeffs(Y^2) == coeffs(Yc^2) == coeffs(Y*Y)

            @test coeffs(Z) == StarAlgebras.mul!(coeffs(W), coeffs(X), coeffs(Y), RG.mstructure)
            @test coeffs(Z) == coeffs(W)
            @test coeffs(Z) == StarAlgebras.mul!(coeffs(W), coeffs(X), coeffs(Y), RGc.mstructure)
            @test coeffs(Z) == coeffs(W)

            StarAlgebras.zero!(W)
            StarAlgebras.fmac!(coeffs(W), coeffs(X), coeffs(Y), RG.mstructure)

            @test coeffs(2*X*Y) == coeffs(StarAlgebras.mul!(W, W, 2))
        end
    end
end


@testset "Group Algebra caching" begin
    A = [:a, :b, :c]
    b = StarAlgebras.Basis{UInt8}(words(A, radius=4))
    k = findfirst(w->length(w)==3, b)-1

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

    @test_throws StarAlgebras.ProductNotDefined StarAlgebras._check(RG.mstructure)

    @test D * D isa AlgebraElement

    @test StarAlgebras._check(RG.mstructure)

    @test all(!iszero, RG.mstructure.table)
end
