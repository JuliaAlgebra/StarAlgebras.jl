# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# Example implementation of Bivariate polynomials
# See MultivariatePolynomials.jl for general implementation
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

function bivariate_algebra()
    exps = ExponentsIterator()
    basis = SA.MappedBasis(exps, monomial, exponents)
    mstr = SA.DiracMStructure(basis, *)
    object = Monomial((0, 0))
    return SA.StarAlgebra(object, mstr)
end
