@testset "FixedBasis caching && allocations" begin
    alph = [:a, :b, :c]
    A★ = FreeWords(alph)
    B = SA.DiracBasis{UInt16}(A★)

    fB = SA.FixedBasis(B; n = nwords(A★, 8), mt = UInt32(nwords(A★, 4)))
    fRG = StarAlgebra(A★, fB)

    k = size(SA.mstructure(basis(fRG)), 1)

    y = spzeros(length(basis(fRG)))
    y[1:k] .= 1
    Y = AlgebraElement(y, fRG)
    @test_broken Y = sum(fRG(basis(fRG)[i]) for i in 1:k)

    @test Y isa AlgebraElement

    @static if v"1.10" ≤ VERSION < v"1.11"
        star(Y)
        star(Y)
        @test (@allocations star(Y)) ≤ 4
    end

    @test supp(Y) == basis(fRG)[1:k]

    @test Y * one(fRG) == Y
    @test one(fRG) * Y == Y

    @test_throws SA.UndefRefError all(!iszero, SA.mstructure(fRG).table)

    @static if v"1.10" ≤ VERSION < v"1.11"
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

    @static if v"1.10" ≤ VERSION < v"1.11"
        YY = deepcopy(Y)

        MA.operate_to!(YY, *, Y, Y)
        @test @allocations(MA.operate_to!(YY, *, Y, Y)) ≤ 25
        @test YY == Y * Y

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
    end
end
