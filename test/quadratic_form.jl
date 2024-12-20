struct Gram{T,B}
    matrix::T
    basis::B
end

function Base.eltype(g::Gram)
    return promote_type(eltype(g.matrix), eltype(eltype(basis(g))))
end
SA.basis(g::Gram) = g.basis
Base.getindex(g::Gram, i, j) = g.matrix[i, j]

@testset "QuadraticForm" begin
    A = let alph = [:a, :b, :c]
        fw = FreeWords(alph)
        SA.StarAlgebra(fw, SA.DiracBasis(fw))
    end

    gbasis = let (id, a, b, c) = A.(Iterators.take(SA.object(A), 4))
        # basis has to be star-invariant:
        bas = 1.0 * [one(A), (a + b) / 2, (a + c) / 2, (b + c) / 2]
        SA.FixedBasis(bas, SA.DiracMStructure(*))
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

end
