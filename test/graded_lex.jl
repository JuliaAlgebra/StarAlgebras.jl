# Example with Graded Lex Ordering
using Test
import StarAlgebras as SA

struct ExponentsIterator end
Base.eltype(::Type{ExponentsIterator}) = NTuple{2,Int}
Base.IteratorSize(::Type{ExponentsIterator}) = Base.IsInfinite()

function Base.iterate(::ExponentsIterator)
    z = (0, 0)
    return z, (z, 0)
end

function Base.iterate(::ExponentsIterator, state)
    z, deg = state
    if iszero(z[2])
        deg += 1
        z = (0, deg)
    else
        z = (z[1] + 1, z[2] - 1)
    end
    return z, (z, deg)
end

function grlex(a::NTuple{2,Int}, b::NTuple{2,Int})
    return isless((sum(a), a), (sum(b), b))
end

struct Monomial
    exponents::NTuple{2,Int}
end

Base.one(::Monomial) = Monomial((0, 0))
Base.:*(a::Monomial, b::Monomial) = Monomial(a.exponents .+ b.exponents)

monomial(exp) = Monomial(exp)
exponents(mono::Monomial) = mono.exponents

SA.comparable(::ExponentsIterator) = grlex

@testset "Graded Lex" begin
    exps = ExponentsIterator()
    basis = SA.MappedBasis(exps, monomial, exponents)
    mstr = SA.DiracMStructure(basis, *)
    object = Monomial((0, 0))
    alg = SA.StarAlgebra(object, mstr)
    @test isone(one(alg))
    a = SA.AlgebraElement(
        SA.SparseCoefficients(
            collect(Iterators.take(exps, 3)),
            [2, -1, 3],
        ),
        alg,
    )
    c = a * a
    @test c.coeffs.values == [4, -4, 12, 1, -6, 9]
    @test c.coeffs.basis_elements == [(0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0)]
end