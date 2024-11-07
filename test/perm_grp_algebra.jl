# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

using Test
using PermutationGroups
import StarAlgebras as SA
@testset "POC: group algebra" begin
    G = PermGroup(perm"(1,2,3,4,5,6)", perm"(1,2)")
    g = Permutation(perm"(1,4,3,6)(2,5)", G)
    h = Permutation(perm"(2,4,5,1)", G)

    db = SA.DiracBasis(G)
    mstr = SA.DiracMStructure(db, *)
    @test mstr(g, h) == SA.SparseCoefficients((g * h,), (1,))

    @test db[g] == g
    @test db[db[g]] == g

    xcfs = SA.SparseCoefficients([one(G), g], [1, -1])
    ycfs = SA.SparseCoefficients([one(G), inv(g)], [1, -1])
    xycfs = SA.SparseCoefficients([one(G), g, inv(g)], [2, -1, -1])

    zcfs = SA.SparseCoefficients([one(G), h], [1, -1])
    xzcfs = SA.SparseCoefficients([one(G), g, h, g * h], [1, -1, -1, 1])

    RG = SA.StarAlgebra(G, mstr)

    x = SA.AlgebraElement(xcfs, RG)
    y = SA.AlgebraElement(ycfs, RG)
    xy = SA.AlgebraElement(xycfs, RG)
    @test x != y
    @test x' == y
    @test x * y == xy

    z = SA.AlgebraElement(zcfs, RG)
    xz = SA.AlgebraElement(xzcfs, RG)
    @test x * z == xz

    # FIXME Broken
#    @testset "Augmented basis" begin
#        ad = SA.AugmentedBasis(db)
#        @test SA.mstructure(ad) == SA.AugmentedMStructure(SA.mstructure(db))
#        @test ad[SA.Augmented(h)] isa SA.Augmented
#        @test sprint(show, ad[SA.Augmented(h)]) == "(-1·()+1·(1,2,4,5))"
#
#        @test !(h in ad)
#        @test SA.Augmented(h) in ad
#
#        IG = SA.StarAlgebra(G, ad)
#
#        axcfs = SA.coeffs(x, basis(IG))
#        aycfs = SA.coeffs(y, basis(IG))
#        azcfs = SA.coeffs(z, basis(IG))
#        ax = SA.AlgebraElement(axcfs, IG)
#        ay = SA.AlgebraElement(aycfs, IG)
#        az = SA.AlgebraElement(azcfs, IG)
#
#        @test coeffs(ax * ay) == SA.coeffs(x * y, basis(IG))
#        @test coeffs(ax * az) == SA.coeffs(x * z, basis(IG))
#        @test SA.aug(ax) == 0
#        @test star(ax) * star(ay) == star(ay) * star(ax)
#
#        @test length(ad) == length(db) - 1
#        @test Set(ad) == Set(SA.Augmented(g) for g in db if !isone(g))
#    end

    @testset "Random elements" begin
        rcfs = SA.SparseCoefficients(rand(G, 10), rand(-2:2, 10))
        r = SA.AlgebraElement(rcfs, RG)
        scfs = SA.SparseCoefficients(rand(G, 10), rand(-2:2, 10))
        s = SA.AlgebraElement(scfs, RG)

        @test SA.aug(r * s) == SA.aug(r) * SA.aug(s)
    end
    @testset "Fixed Basis" begin
        m = PermutationGroups.order(UInt16, G)
        fb = SA.FixedBasis{eltype(G),typeof(m)}(collect(G))

        @test fb[fb[g]] == g

        fRG = SA.StarAlgebra(G, SA.MTable(fb, (m, m)))

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

        a, b = let mt = SA.mstructure(fRG).table
            count(i -> isassigned(mt, i), eachindex(mt)), length(mt)
        end
        @test a ≤ b
    end

    @testset "FiniteSupportBasis" begin
        S1 = unique!(rand(G, 7))
        S = unique!([S1; [a * b for a in S1 for b in S1]])
        subb = SA.SubBasis(S, db)
        smstr = SA.mstructure(subb)
        @test smstr(1, 2).basis_elements[1] == subb[1] * subb[2]

        sbRG = SA.StarAlgebra(G, subb)

        x = let z = zeros(Int, length(SA.basis(sbRG)))
            z[1:length(S1)] .= rand(-1:1, length(S1))
            SA.AlgebraElement(z, sbRG)
        end

        y = let z = zeros(Int, length(SA.basis(sbRG)))
            z[1:length(S1)] .= rand(-1:1, length(S1))
            SA.AlgebraElement(z, sbRG)
        end

        dx = SA.AlgebraElement(SA.coeffs(x, SA.basis(RG)), RG)
        dy = SA.AlgebraElement(SA.coeffs(y, SA.basis(RG)), RG)

        @test dx + dy == SA.AlgebraElement(SA.coeffs(x + y, SA.basis(RG)), RG)

        @test dx * dy == SA.AlgebraElement(SA.coeffs(x * y, SA.basis(RG)), RG)

        a = SA.AlgebraElement([2], SA.StarAlgebra(G, SA.SubBasis([g], db)))
        b = SA.AlgebraElement([-3], SA.StarAlgebra(G, SA.SubBasis([h], db)))
        # `Base.+` assumes that using the basis of the first argument will suffice
        # We should redefine `Base.:+(a::SubBasis, b::SubBasis)` to first
        # convert `a` and `b` to their implicit basis equivalent and then
        # do `+` and then convert the result back
        # `MultivariateBases` defines an `implicit` function.
        # Why not having an `explicit` as well ?
        # My dream implementation would be
        # Base.:+(a::SubBasis, b::SubBasis) = explicit(implicit(a) + implicit(b))
        # so we just need to implement `implicit` and `explicit` 👼
        @test_broken SA.explicit(SA.implicit(a)) == a
        @test_broken SA.explicit(SA.implicit(b)) == b
        @test_broken a + b == SA.explicit(SA.implicit(a) + SA.implicit(b))
        @test_broken a * b == SA.explicit(SA.implicit(a) * SA.implicit(b))
    end
end
