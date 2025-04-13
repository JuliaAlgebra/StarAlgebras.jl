# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

function _alloc_test(output, op, a, b, allocs)
    MA.operate_to!(output, op, a, b)
    expected = op(a, b)
    @test @allocations(MA.operate_to!(output, op, a, b)) <= allocs
    @test output == expected
end

function _test_op(op, a, b)
    result = @inferred op(a, b)
    @test typeof(result) == MA.promote_operation(op, typeof(a), typeof(b))
    return result
end

@testset "FixedBasis caching && allocations" begin
    alph = [:a, :b, :c]
    A★ = FreeWords(alph)
    B = SA.DiracBasis(A★)

    fB = SA.FixedBasis(B; n = nwords(A★, 8), mt = UInt32(nwords(A★, 4)))
    fRG = StarAlgebra(A★, fB)

    k = size(SA.mstructure(basis(fRG)), 1)

    y = spzeros(length(basis(fRG)))
    y[1:k] .= 1
    Y = AlgebraElement(y, fRG)
    @test Y == sum(fRG(basis(fRG)[i]) for i in 1:k)

    @test Y isa AlgebraElement

    @static if v"1.10" ≤ VERSION
        star(Y)
        star(Y)
        @test (@allocations star(Y)) ≤ 6
    end

    @test SA.supp(Y) == basis(fRG)[1:k]

    @test Y * one(fRG) == Y
    @test one(fRG) * Y == Y

    @test_throws SA.UndefRefError all(!iszero, SA.mstructure(fRG).table)
    mstr = deepcopy(SA.mstructure(fRG))
    SA.complete!(mstr)
    @test all(!iszero(mstr.table))
    @test_throws SA.UndefRefError all(!iszero, SA.mstructure(fRG).table)

    @static if v"1.10" ≤ VERSION
        @test (@allocations Y * Y) > 2(k^2 - 2 * k)
        @test Y * Y isa AlgebraElement
        res = SA._preallocate_output(*, Y, Y)
        @test @allocations(SA._preallocate_output(*, Y, Y)) ≤ 10
        MA.operate_to!(res, *, Y, Y) # makes sure res is of right size
        @test (@allocations MA.operate_to!(res, *, Y, Y)) ≤ 1
    end

    @test all(!iszero, SA.mstructure(fRG).table)

    @static if v"1.10" ≤ VERSION
        YY = deepcopy(Y)
        _alloc_test(YY, *, Y, Y, 25)
        _alloc_test(YY, +, Y, Y, 0)
        _alloc_test(YY, -, Y, Y, 0)

        # SparseArrays calls `Base.unalias` which allocates:
        _alloc_test(YY, +, YY, Y, 4)
        _alloc_test(YY, +, Y, YY, 4)
    end
end

@testset "tuple" begin
    alph = [:a, :b, :c]
    A★ = FreeWords(alph)
    B = SA.DiracBasis(A★)

    fB = SA.FixedBasis(B; n = nwords(A★, 2), mt = UInt32(nwords(A★, 2)))
    fRG = StarAlgebra(A★, fB)

    y = SA.SparseCoefficients(
        [first(fB)],
        [1],
    )
    Y = AlgebraElement(y, fRG)

    z = SA.SparseCoefficients(
        (first(fB),),
        (1,),
    )
    Z = AlgebraElement(z, fRG)
    @test _test_op(+, Z, Z) == _test_op(*, 2, Z)
    @test _test_op(+, Z, Z) == _test_op(+, Y, Y)
    @test _test_op(+, Y, Z) == _test_op(+, Y, Y)
    @test _test_op(+, Z, Y) == _test_op(+, Y, Y)
    @test _test_op(-, Z, Z) == _test_op(*, 0, Z)
    @test _test_op(-, Z, Z) == _test_op(-, Y, Z)

    for X in [Y, Z]
        c = coeffs(X)
        res = 2 .* c
        @test c .* 2 == res
        @test c .+ 1 == res
        @test 1 .+ c == res
        @test MA.Zero() .+ c == c
        @test c .+ MA.Zero() == c
        err = ArgumentError("Cannot broadcast `StarAlgebras.SparseCoefficients` with another array of different type")
        @test_throws err c .+ ones(3)
    end
end
