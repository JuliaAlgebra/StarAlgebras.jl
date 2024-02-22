using PermutationGroups
using Random

SA.star(x::Number) = x'
SA.star(g::PermutationGroups.AP.AbstractPermutation) = inv(g)

@testset "POC: group algebra" begin
    G = PermGroup(perm"(1,2,3,4,5,6)", perm"(1,2)")
    g = Permutation(perm"(1,4,3,6)(2,5)", G)
    h = Permutation(perm"(2,4,5,1)", G)

    db = SA.DiracBasis{UInt32}(G)
    @test SA.mstructure(db) == SA.DiracMStructure(*)
    @test SA.mstructure(db)(g, h) == SA.Dirac(g * h)

    @test db[g] isa SA.Dirac
    @test_throws MethodError db[db[g]]

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
        @test ad[SA.AugmentedDirac(h)] isa SA.AugmentedDirac

        IG = SA.StarAlgebra(G, ad)

        axcfs = SA.coeffs(x, basis(IG))
        aycfs = SA.coeffs(y, basis(IG))
        azcfs = SA.coeffs(z, basis(IG))
        ax = SA.AlgebraElement(axcfs, IG)
        ay = SA.AlgebraElement(aycfs, IG)
        az = SA.AlgebraElement(azcfs, IG)

        @test coeffs(ax * ay) == SA.coeffs(x * y, basis(IG))
        @test coeffs(ax * az) == SA.coeffs(x * z, basis(IG))
    end

    @testset "Random elements" begin
        rcfs = SA.SparseCoefficients(rand(G, 10), rand(-2:2, 10))
        r = SA.AlgebraElement(rcfs, RG)
        scfs = SA.SparseCoefficients(rand(G, 10), rand(-2:2, 10))
        s = SA.AlgebraElement(scfs, RG)

        @test aug(r * s) == aug(r) * aug(s)
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

        @test aug(fr) == aug(r)
        @test aug(fs) == aug(s)
        @test aug(fr * fs) == aug(fr) * aug(fs)

        rs_cfs = coeffs(r * s, basis(fRG))
        @test rs_cfs isa AbstractVector
        @test fr * fs == SA.AlgebraElement(rs_cfs, fRG)

        a, b = let mt = SA.mstructure(basis(fRG)).table
            count(i -> isassigned(mt, i), eachindex(mt)), length(mt)
        end
        @test a ≤ b
    end

end
