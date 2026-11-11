# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

using Test
using PermutationGroups
import Random
import StarAlgebras as SA
import MutableArithmetics as MA
if !@isdefined aug
    include("group_utils.jl")
end

@testset "POC: group algebra" begin
    G = PermGroup(perm"(1,2,3,4,5,6)", perm"(1,2)")
    g = Permutation(perm"(1,4,3,6)(2,5)", G)
    h = Permutation(perm"(2,4,5,1)", G)

    db = SA.DiracBasis(G)
    mstr = SA.DiracMStructure(db, *)
    @test mstr(g, h) == SA.SparseCoefficients((g * h,), (1,))

    @test db[g] == g
    @test db[db[g]] == g

    xcfs = SA.sparse_coefficients([one(G), g], [1, -1])
    ycfs = SA.sparse_coefficients([one(G), inv(g)], [1, -1])
    xycfs = SA.sparse_coefficients([one(G), g, inv(g)], [2, -1, -1])

    zcfs = SA.sparse_coefficients([one(G), h], [1, -1])
    xzcfs = SA.sparse_coefficients([one(G), g, h, g * h], [1, -1, -1, 1])

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

    @testset "Augmented basis" begin
        ad = AugmentedBasis(db)
        @test SA.mstructure(ad) == AugmentedMStructure(ad)
        @test ad[Augmented(h)] isa Augmented
        @test sprint(show, ad[Augmented(h)]) == "(-1·()+1·(1,2,4,5))"

        @test !(h in ad)
        @test Augmented(h) in ad

        IG = SA.StarAlgebra(G, ad)

        axcfs = SA.coeffs(x, basis(IG))
        aycfs = SA.coeffs(y, basis(IG))
        azcfs = SA.coeffs(z, basis(IG))
        ax = SA.AlgebraElement(axcfs, IG)
        ay = SA.AlgebraElement(aycfs, IG)
        az = SA.AlgebraElement(azcfs, IG)

        @test coeffs(ax * ay) == SA.coeffs(x * y, basis(IG))
        @test coeffs(ax * az) == SA.coeffs(x * z, basis(IG))
        @test aug(ax) == 0
        @test star(ax) * star(ay) == star(ay) * star(ax)

        @test length(ad) == length(db) - 1
        @test Set(ad) == Set(Augmented(g) for g in db if !isone(g))
    end

    @testset "Random elements (seed=$seed)" for seed in 0:20
        Random.seed!(seed)
        rcfs = SA.sparse_coefficients(rand(G, 10), rand(-2:2, 10))
        r = SA.AlgebraElement(rcfs, RG)
        scfs = SA.sparse_coefficients(rand(G, 10), rand(-2:2, 10))
        s = SA.AlgebraElement(scfs, RG)

        @test aug(r * s) == aug(r) * aug(s)
    end
    @testset "Fixed Basis" begin
        m = PermutationGroups.order(UInt16, G)
        fb = SA.FixedBasis{eltype(G),typeof(m)}(collect(G))

        @test fb[fb[g]] == g

        fRG = SA.StarAlgebra(G, SA.MTable(fb, (m, m)))

        @testset "seed=$seed" for seed in 0:20
            Random.seed!(seed)
            rcfs = SA.sparse_coefficients(rand(G, 10), rand(-2:2, 10))
            r = SA.algebra_element(rcfs, RG)
            scfs = SA.sparse_coefficients(rand(G, 10), rand(-2:2, 10))
            s = SA.algebra_element(scfs, RG)

            @test coeffs(r, basis(fRG)) isa SparseVector

            fr = SA.AlgebraElement(coeffs(r, basis(fRG)), fRG)
            fs = SA.AlgebraElement(coeffs(s, basis(fRG)), fRG)

            @test aug(fr) == aug(r)
            @test aug(fs) == aug(s)
            @test aug(fr * fs) == aug(fr) * aug(fs)

            @test coeffs(r * s, basis(fRG)) isa AbstractVector
            @test fr * fs == SA.AlgebraElement(coeffs(r * s, basis(fRG)), fRG)
        end

        a, b = let mt = SA.mstructure(fRG).table
            count(i -> isassigned(mt, i), eachindex(mt)), length(mt)
        end
        @test a ≤ b

        @test Random.seed!(0) isa Any
        rcfs = SA.sparse_coefficients(rand(G, 10), rand(-2:2, 10))
        r = SA.AlgebraElement(rcfs, RG)
        fr = SA.AlgebraElement(coeffs(r, basis(fRG)), fRG)

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
        @test only(smstr(1, 2).basis_elements) == subb[subb[1]*subb[2]]
        @test only(smstr(1, 2, eltype(subb)).basis_elements) ==
              subb[1] * subb[2]

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

        @test SA.implicit_basis(SA.basis(sbRG)) === SA.basis(RG)
        @test SA.implicit_basis(SA.basis(RG)) === SA.basis(RG)
        @test SA.implicit(x) == dx
        @test parent(SA.implicit(x)) == RG
        @test SA.implicit(dx) === dx
        @test typeof(SA.implicit(x)) ==
              MA.promote_operation(SA.implicit, typeof(x))
        @test MA.promote_operation(SA.implicit_basis, typeof(SA.basis(sbRG))) ==
              typeof(SA.basis(RG))
        @test MA.promote_operation(SA.implicit_basis, typeof(SA.basis(RG))) ==
              typeof(SA.basis(RG))
        ts = SA.Term(sbRG, 1, 2)
        @test SA.implicit(ts).algebra == RG
        @test SA.basis_element(SA.implicit(ts)) == SA.basis_element(ts)
        @test SA.coefficient(SA.implicit(ts)) == 2
        @test SA.implicit(SA.implicit(ts)) === SA.implicit(ts)

        # Arithmetic between a term of the algebra of the parent basis and an
        # element of the algebra of a `SubBasis`
        t = SA.Term(RG, SA.basis(RG)[S1[1]], 2)
        dt = 2 * RG(S1[1])
        for op in (+, -, *)
            @test op(t, x) == op(dt, dx)
            @test op(x, t) == op(dx, dt)
            @test parent(op(t, x)) == RG
        end

        # The constant term may not be in `subb` so the sum is in the algebra of
        # the parent basis
        for op in (+, -)
            @test op(x, 2) == op(dx, 2)
            @test op(2, x) == op(2, dx)
            @test parent(op(x, 2)) == RG
            @test typeof(op(x, 2)) == MA.promote_operation(op, typeof(x), Int)
        end

        # The parents cannot be promoted to a common parent
        RG2 = SA.StarAlgebra(G, SA.DiracMStructure(db, (a, b) -> a * b))
        @test_throws ArgumentError MA.promote_operation(
            *,
            typeof(dx),
            typeof(SA.AlgebraElement(SA.coeffs(dx), RG2)),
        )

        # Arithmetic between elements of the algebra of a `SubBasis` and
        # elements of the algebra of its parent basis
        @test x == dx
        @test hash(x) == hash(dx)
        @test x != dy
        for op in (+, -, *)
            @test op(x, dy) == op(dx, dy)
            @test parent(op(x, dy)) == RG
            @test op(dx, y) == op(dx, dy)
            @test parent(op(dx, y)) == RG
            @test typeof(op(x, dy)) ==
                  MA.promote_operation(op, typeof(x), typeof(dy))
            @test typeof(op(dx, y)) ==
                  MA.promote_operation(op, typeof(dx), typeof(y))
        end

        # The products of the first `length(S1)` elements of `subb` are in
        # `subb` so they are given by an `MTable` and the arithmetic stays in
        # the algebra of the `SubBasis`
        tRG = SA.StarAlgebra(G, SA.MTable(smstr, (length(S1), length(S1))))
        tx = SA.AlgebraElement(SA.coeffs(x), tRG)
        ty = SA.AlgebraElement(SA.coeffs(y), tRG)
        @test parent(SA.implicit(tx)) == RG
        @test SA.implicit(tx) == dx
        @test SA.promote_operand(tx) === tx
        @test SA.promote_operand(x) == dx
        @test parent(SA.promote_operand(x)) == RG
        @test SA.promote_operand(dx) === dx
        for op in (+, -, *)
            @test parent(op(tx, ty)) === tRG
            @test typeof(op(tx, ty)) ==
                  MA.promote_operation(op, typeof(tx), typeof(ty))
            @test SA.AlgebraElement(SA.coeffs(op(tx, ty), SA.basis(RG)), RG) ==
                  op(dx, dy)
        end

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
            # Elements of the algebras of two different `SubBasis` of the
            # same parent basis
            z = SA.AlgebraElement(SA.coeffs(dx, SA.basis(sbRG)), sbRG)
            @test parent(x + z) == RG
            @test x + z ==
                  SA.AlgebraElement(SA.coeffs(x, SA.basis(RG)), RG) + dx
            @test eltype(one(Float64, sRG)) == Float64
            @test isone(one(x))
            @test coeffs(one(x)) == coeffs(one(sRG))
        end
    end
    @testset "Algebra Elements Basis (seed = $seed)" for seed in 0:200
        Random.seed!(seed)
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

        ar = SA.AlgebraElement(
            coeffs(rcfs, SA.DiracBasis(elts), basis(aRG)),
            aRG,
        )
        as = SA.AlgebraElement(
            coeffs(scfs, SA.DiracBasis(elts), basis(aRG)),
            aRG,
        )

        @test aug(ar) == 2*ar(first(ab)) # one(RG)
        @test aug(as) == 2as(first(ab)) # one(RG)
        @test aug(ar + as) == aug(ar) + aug(as)

        @test coeffs(ar + as, basis(aRG)) isa AbstractVector

        @test_throws ArgumentError one(aRG)
        if iszero(ar)
            @test !isone(ar)
        else
            @test_throws ArgumentError isone(ar)
        end
    end
end
