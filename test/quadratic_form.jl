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

struct SubBasis{T,I,V<:AbstractVector{I},B<:SA.ImplicitBasis{T,I}} <: SA.ExplicitBasis{Float64,Int}
    implicit::B
    indices::V
end
Base.length(b::SubBasis) = length(b.indices)
function Base.iterate(b::SubBasis, args...)
    elem_state = iterate(b.indices, args...)
    if isnothing(elem_state)
        return
    end
    return b.implicit[elem_state[1]], elem_state[2]
end

@testset "Int -> Float basis" begin
    implicit = SA.MappedBasis(NaturalNumbers(), float, Int)
    explicit = SubBasis(implicit, 1:3)
    m = Bool[
        true  false true
        false true  false
        true  false true
    ]
    Q = SA.QuadraticForm(Gram(m, explicit))
    A = SA.StarAlgebra(nothing, implicit)
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
