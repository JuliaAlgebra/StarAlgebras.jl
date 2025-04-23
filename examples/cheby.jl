# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# This file tests the implementation of a Chebyshev basis using MappedBasis.

# Define a struct to represent a Chebyshev polynomials
struct ChebyPoly
    n::Int
end

function Base.isless(a::ChebyPoly, b::ChebyPoly)
    return a.n < b.n
end

SA.star(a::ChebyPoly) = a

struct ChebyMStruct{T,I,B<:SA.AbstractBasis{T,I}} <: SA.MultiplicativeStructure{T,I}
    basis::B
end

function (m::ChebyMStruct)(a::ChebyPoly, b::ChebyPoly, ::Type{U}) where {U}
    return m(a.n, b.n, U)
end

function (m::ChebyMStruct)(a::Int, b::Int, ::Type{ChebyPoly})
    return SA.map_keys(Base.Fix1(getindex, m), m(a, b, Int))
end

function (m::ChebyMStruct)(a::Int, b::Int, ::Type{Int})
    return SA.SparseCoefficients(
        (abs(a - b), a + b),
        (1 // 2, 1 // 2),
    )
end

# Define the mapping from `Int` to `ChebyPoly`
function map_to_chebyshev(n::Int)
    return ChebyPoly(n)
end

# Define the inverse mapping from `ChebyPoly`s to `NaturalNumber`
function map_from_chebyshev(p::ChebyPoly)
    return p.n
end

# Create a `MappedBasis` that maps `NaturalNumbers` to `ChebyPoly`s
function cheby_basis()
    return SA.MappedBasis(NaturalNumbers(), map_to_chebyshev, map_from_chebyshev)
end