struct Gram{U,T,B}
    matrix::T
    basis::B
end

function Gram{U}(matrix, basis) where {U}
    return Gram{U,typeof(matrix),typeof(basis)}(matrix, basis)
end

function Gram(matrix, basis)
    U = promote_type(eltype(matrix), eltype(eltype(basis)))
    return Gram{U}(matrix, basis)
end

Base.eltype(::Gram{U}) where {U} = U

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

struct IntToFloat <: SA.ImplicitBasis{Float64,Int} end
SA.mstructure(::IntToFloat) = SA.DiracMStructure(*)
Base.first(::IntToFloat) = 1.0
Base.getindex(::IntToFloat, i::Int) = convert(Float64, i)
Base.getindex(::IntToFloat, i::Float64) = convert(Int, i)

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

@testset "Dummy" begin
    implicit = IntToFloat()
    explicit = SubBasis(implicit, 1:3)
    m = Bool[
        true  false true
        false true  false
        true  false true
    ]
    Q = SA.QuadraticForm(Gram{Int}(m, explicit))
    A = SA.StarAlgebra(nothing, implicit)
    @test A(Q) == SA.AlgebraElement(
        SA.SparseCoefficients(
            [1.0, 3.0, 4.0, 9.0],
            [1, 2, 1, 1],
        ),
        A,
    )
    mt = SA.MTable(float.(1:6), SA.DiracMStructure(*), (0, 0))
    @test mt(2.0, 3.0) == SA.SparseCoefficients([6.0], [1])
    @test mt(2.0, 3.0) == mt(2, 3)
end
