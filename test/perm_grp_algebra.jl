# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

using Test
using PermutationGroups
import Random
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

    @test isone(one(RG))
    @test x * one(RG) == x
    @test one(RG) * x == x
    @test eltype(one(Float64, RG)) == Float64
    @test isone(one(x))
    @test coeffs(one(x)) == coeffs(one(RG))

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

        @test isone(one(fRG))
        @test fr * one(fRG) == fr
        @test one(fRG) * fr == fr
        @test eltype(one(Float64, fRG)) == Float64
        @test isone(one(fr))
        @test coeffs(one(fr)) == coeffs(one(fRG))
    end

    @testset "SubBasis" begin
        # If we're unlucky, `a^3` might belong to the basis.
        # We fix the seed to be sure we are never in that case.
        Random.seed!(0)
        S1 = unique!(rand(G, 7))
        S = unique!([S1; [a * b for a in S1 for b in S1]])
        subb = SA.sub_basis(db, S)
        a = S1[1]
        @test subb[a] == 1
        @test a in subb
        @test isnothing(get(subb, a^3, nothing))
        @test_throws KeyError(a^3) subb[a^3]
        @test !(a^3 in subb)
        @test collect(subb) == S
        smstr = SA.DiracMStructure(subb, *)
        @test only(smstr(1, 2).basis_elements) == subb[subb[1] * subb[2]]
        @test only(smstr(1, 2, eltype(subb)).basis_elements) == subb[1] * subb[2]

        sbRG = SA.StarAlgebra(G, smstr)

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

        if !(one(G) in subb)
            @test_throws ArgumentError one(sbRG)
        end

        S2 = unique([S; one(G)])
        subb2 = SA.sub_basis(db, S2)
        let sRG = SA.StarAlgebra(G, subb2)
            x = let z = spzeros(Int, length(SA.basis(sRG)))
                z[rand(1:length(S2), 10)] += rand(-1:1, 10)
                SA.AlgebraElement(z, sRG)
            end
            @test isone(one(sRG))
            @test x * one(sRG) == x
            @test one(sRG) * x == x
            @test eltype(one(Float64, sRG)) == Float64
            @test isone(one(x))
            @test coeffs(one(x)) == coeffs(one(sRG))
        end
    end
    @testset "Algebra Elements Basis" begin
        S1 = unique!(collect(Iterators.take(G, 10)))
        S = unique!([a * b for a in S1 for b in S1])

        elts = [RG(g) - RG(1) for g in S]
        elts[begin] = RG(S[2]) + RG(1)
        ab = SA.FixedBasis(elts)
        # basis looks like this: [s₁ + e, s₁ - e, s₂ - e, s₃ - e, …]
        # in particular it doesn't contain 1!

        @test ab[ab[elts[2]]] == elts[2]

        rcfs = SA.SparseCoefficients(rand(elts[1:length(S1)], 5), rand(-2:2, 5))
        scfs = SA.SparseCoefficients(rand(elts[1:length(S1)], 5), rand(-2:2, 5))

        @test coeffs(rcfs, SA.DiracBasis(elts), ab) isa SparseVector

        aRG = SA.StarAlgebra(RG, SA.MTable(ab, (length(S1), length(S1))))

        ar = SA.AlgebraElement(coeffs(rcfs, SA.DiracBasis(elts), basis(aRG)), aRG)
        as = SA.AlgebraElement(coeffs(scfs, SA.DiracBasis(elts), basis(aRG)), aRG)

        @test SA.aug(ar) == 2*ar(first(ab)) # one(RG)
        @test SA.aug(as) == 2as(first(ab)) # one(RG)
        @test SA.aug(ar + as) == SA.aug(ar) + SA.aug(as)

        @test coeffs(ar + as, basis(aRG)) isa AbstractVector

        @test_throws ArgumentError one(aRG)
        @test_throws ArgumentError isone(ar)
    end
end
