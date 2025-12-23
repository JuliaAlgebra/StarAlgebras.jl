# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# Example implementation of Univariate polynomials
# See MultivariatePolynomials.jl for general implementation
struct UnivariateIterator end
Base.eltype(::Type{UnivariateIterator}) = Int
Base.IteratorSize(::Type{UnivariateIterator}) = Base.IsInfinite()

function Base.iterate(::UnivariateIterator)
    z = 0
    return z, z
end

function Base.iterate(::UnivariateIterator, z)
    z += 1
    return z, z
end

struct UniMono
    first_variable::Bool
    exponent::Int
end

Base.one(m::UniMono) = UniMono(m.first_variable, 0)
function Base.:*(a::UniMono, b::UniMono)
    @assert a.first_variable == b.first_variable
    UniMono(a.first_variable, a.exponent + b.exponent)
end

unimono(first, exp) = UniMono(first, exp)
exponent(mono::UniMono) = mono.exponent

SA.comparable(::UnivariateIterator) = isless

function univariate_algebra(first_variable::Bool)
    exps = UnivariateIterator()
    basis = SA.MappedBasis(exps, Base.Fix1(unimono, first_variable), exponent)
    mstr = SA.DiracMStructure(basis, *)
    object = UniMono(first_variable, 0)
    return SA.StarAlgebra(object, mstr)
end

function bivariate_lift(first, deg)
    return first ? (deg, 0) : (0, deg)
end

UniType = typeof(SA.MappedBasis(UnivariateIterator(), Base.Fix1(unimono, true), exponent))
function SA.promote_basis_with_maps(a::UniType, b::UniType)
    if a.map.x == b.map.x
        return (a, nothing), (b, nothing)
    end
    bi_alg = bivariate_algebra()
    bi_bas = SA.basis(bi_alg)
    _a = (bi_bas, Base.Fix1(bivariate_lift, a.map.x))
    _b = (bi_bas, Base.Fix1(bivariate_lift, b.map.x))
    return _a, _b
end

function SA.promote_object(m::UniMono, mstr, map::Base.Fix1)
    return Monomial(map(m.exponent))
end
