@testset "POC: group algebra" begin
    G = PermGroup(perm"(1,2,3,4,5,6)", perm"(1,2)")
    g = Permutation(perm"(1,4,3,6)(2,5)", G)
    h = Permutation(perm"(2,4,5,1)", G)

    db = SA.DiracBasis{UInt32}(G)
    @test SA.mstructure(db) == SA.DiracMStructure(*)
    @test SA.mstructure(db)(g, h) == SA.SparseCoefficients((g * h,), (1,))

    @test db[g] == g
    @test db[db[g]] == g

    xcfs = SA.SparseCoefficients([one(G), g], [1, -1])
    ycfs = SA.SparseCoefficients([one(G), inv(g)], [1, -1])
    xycfs = SA.SparseCoefficients([one(G), g, inv(g)], [2, -1, -1])

    zcfs = SA.SparseCoefficients([one(G), h], [1, -1])
    xzcfs = SA.SparseCoefficients([one(G), g, h, g * h], [1, -1, -1, 1])

    RG = SA.StarAlgebra(G, db)

    x = SA.AlgebraElement(xcfs, RG)
    y = SA.AlgebraElement(ycfs, RG)
    xy = SA.AlgebraElement(xycfs, RG)
    @test x != y
    @test x' == y
    @test x * y == xy

    z = SA.AlgebraElement(zcfs, RG)
    xz = SA.AlgebraElement(xzcfs, RG)
    @test x * z == xz

    @testset "Augmented basis" begin

        ad = SA.AugmentedBasis(db)
        @test SA.mstructure(ad) == SA.AugmentedMStructure(SA.mstructure(db))
        @test ad[SA.Augmented(h)] isa SA.Augmented

        IG = SA.StarAlgebra(G, ad)

        axcfs = SA.coeffs(x, basis(IG))
        aycfs = SA.coeffs(y, basis(IG))
        azcfs = SA.coeffs(z, basis(IG))
        ax = SA.AlgebraElement(axcfs, IG)
        ay = SA.AlgebraElement(aycfs, IG)
        az = SA.AlgebraElement(azcfs, IG)

        @test coeffs(ax * ay) == SA.coeffs(x * y, basis(IG))
        @test coeffs(ax * az) == SA.coeffs(x * z, basis(IG))
        @test SA.aug(ax) == 0
    end

    @testset "Random elements" begin
        rcfs = SA.SparseCoefficients(rand(G, 10), rand(-2:2, 10))
        r = SA.AlgebraElement(rcfs, RG)
        scfs = SA.SparseCoefficients(rand(G, 10), rand(-2:2, 10))
        s = SA.AlgebraElement(scfs, RG)

        @test SA.aug(r * s) == SA.aug(r) * SA.aug(s)
    end
    @testset "Fixed Basis" begin
        m = PermutationGroups.order(UInt16, G)
        fb = SA.FixedBasis(collect(G), SA.DiracMStructure(*), (m, m))

        @test fb[fb[g]] == g

        fRG = SA.StarAlgebra(G, fb)

        rcfs = SA.SparseCoefficients(rand(G, 10), rand(-2:2, 10))
        r = SA.AlgebraElement(rcfs, RG)
        scfs = SA.SparseCoefficients(rand(G, 10), rand(-2:2, 10))
        s = SA.AlgebraElement(scfs, RG)

        @test coeffs(r, basis(fRG)) isa SparseVector

        fr = SA.AlgebraElement(coeffs(r, basis(fRG)), fRG)
        fs = SA.AlgebraElement(coeffs(s, basis(fRG)), fRG)

        @test SA.aug(fr) == SA.aug(r)
        @test SA.aug(fs) == SA.aug(s)
        @test SA.aug(fr * fs) == SA.aug(fr) * SA.aug(fs)

        @test coeffs(r * s, basis(fRG)) isa AbstractVector
        @test fr * fs == SA.AlgebraElement(coeffs(r * s, basis(fRG)), fRG)

        a, b = let mt = SA.mstructure(basis(fRG)).table
            count(i -> isassigned(mt, i), eachindex(mt)), length(mt)
        end
        @test a ≤ b
    end
end
