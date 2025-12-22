# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

@testset "QuadraticForm" begin
    A = let alph = [:a, :b, :c]
        fw = FreeWords(alph)
        SA.StarAlgebra(fw, SA.DiracBasis(fw))
    end

    gbasis = let (id, a, b, c) = A.(Iterators.take(SA.object(A), 4))
        # basis has to be star-invariant:
        bas = 1.0 * [one(A), (a + b) / 2, (a + c) / 2, (b + c) / 2]
        SA.FixedBasis(bas)
    end

    m = [
        π 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
    ]
    Q = SA.QuadraticForm(Gram(m, gbasis))
    b = basis(Q)
    @test A(Q) == π * b[1] * b[1]

    m = [
        0 1//2 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
    ]
    Q = SA.QuadraticForm(Gram(m, gbasis))
    b = basis(Q)
    @test A(Q) == 1 // 2 * b[1] * b[2]

    m = [
        0 0 0 0
        0 π -1 0
        0 0 0 0
        0 0 0 0
    ]
    Q = SA.QuadraticForm(Gram(m, gbasis))
    b = basis(Q)
    @test A(Q) == π * b[2] * b[2] - b[2] * b[3]

    m = [
        0 0 0 0
        0 0 0 0
        0 π 0 0
        0 0 1 0
    ]
    Q = SA.QuadraticForm(Gram(m, gbasis))
    b = basis(Q)
    @test A(Q) == π * b[3] * b[2] + b[4] * b[3]

    Q = SA.QuadraticForm{SA.star}(Gram(m, gbasis))
    @test A(Q) == π * b[3]' * b[2] + b[4]' * b[3]

    m = ones(Int, 4, 4)
    Q = SA.QuadraticForm(Gram(m, gbasis))
    @test A(Q) == sum(bi * bj for bi in gbasis for bj in gbasis)

    Q = SA.QuadraticForm{SA.star}(Gram(m, gbasis))
    @test A(Q) == sum(bi' * bj for bi in gbasis for bj in gbasis)

    m = [
        0 0 0 0
        0 0 0 0
        0 π 0 0
        0 0 1 0
    ]
    Q = SA.QuadraticForm(Gram(m, gbasis))
    b = basis(Q)
    @test A(Q) == π * b[3] * b[2] + b[4] * b[3]
end

@testset "Int -> Float basis" begin
    limited = SA.MappedBasis(1:2, float, Int)
    @test !(3.0 in limited)
    @test !haskey(limited, 3)
    @test collect(limited) == [1.0, 2.0]
    implicit = SA.MappedBasis(NaturalNumbers(), float, Int)
    @test 3.0 in implicit
    @test haskey(implicit, 3)
    explicit = SA.SubBasis(implicit, 1:3)
    @test 3.0 in explicit
    @test collect(explicit) == [1.0, 2.0, 3.0]
    @test haskey(explicit, 3)
    m = Bool[
        true  false true
        false true  false
        true  false true
    ]
    Q = SA.QuadraticForm(Gram(m, explicit))
    A = SA.StarAlgebra(PlaceholderObject(), implicit)
    @test A(Q) == SA.AlgebraElement(
        SA.SparseCoefficients(
            [1.0, 3.0, 4.0, 9.0],
            [1, 2, 1, 1],
        ),
        A,
    )
    mt = SA.MTable(implicit, (0, 0))
    @test mt(2.0, 3.0) == SA.SparseCoefficients([6.0], [1])
    @test mt(2.0, 3.0) == mt(2, 3)
end

@testset "Chebyshev basis" begin
    implicit = cheby_basis()
    mstr = ChebyMStruct(implicit)
    mt = SA.MTable(mstr, (0, 0))
    sub = SA.SubBasis(implicit, 1:3)
    test_vector_interface(sub)
    fixed = SA.FixedBasis(implicit; n = 3)
    test_vector_interface(fixed)
    a = ChebyPoly(2)
    b = ChebyPoly(3)
    expected = SA.SparseCoefficients((ChebyPoly(1), ChebyPoly(5)), (1 // 2, 1 // 2))
    for mult in [mstr, mt]
        @test mult(a, b) == expected
        @test mult(a, b, Int) == mult(2, 3)
        @test mult(a, b) == mult(2, 3, ChebyPoly)
        m = [
            2      1 // 2 0
            1 // 2 0      2
            0      2      1
        ]
        A = SA.StarAlgebra(PlaceholderObject(), mult)
        for explicit in [sub, fixed]
            Q = SA.QuadraticForm(Gram(m, explicit))
            @test A(Q) == SA.AlgebraElement(
                SA.SparseCoefficients(
                    [0, 1, 2, 3, 5, 6],
                    [3//2, 5//2, 1, 1//2, 2, 1//2],
                ),
                A,
            )
        end
    end
end
