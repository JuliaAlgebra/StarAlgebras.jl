# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

mutable struct FixedBasis{T,I,V<:AbstractVector{T}} <:
               ExplicitBasis{T,I}
    elts::V
    relts::Dict{T,I}
    starof::Vector{I}
end

function FixedBasis{T,I}(elts::AbstractVector{T}) where {T,I}
    relts = Dict(b => I(idx) for (idx, b) in pairs(elts))
    starof = [relts[star(x)] for x in elts]
    return FixedBasis{T,I,typeof(elts)}(elts, relts, starof)
end

FixedBasis(elts::AbstractVector{T}) where {T} = FixedBasis{T,keytype(elts)}(elts)

function FixedBasis{T,I}(basis::AbstractBasis{T}; n::Integer) where {T,I}
    return FixedBasis{T,I}(collect(Iterators.take(basis, n)))
end

FixedBasis(basis::AbstractBasis{T}; n::Integer) where {T} = FixedBasis{T,typeof(n)}(basis; n)

Base.in(x, b::FixedBasis) = haskey(b.relts, x)
Base.getindex(b::FixedBasis{T}, x::T) where {T} = b.relts[x]
Base.getindex(b::FixedBasis, i::Integer) = b.elts[i]

Base.IteratorSize(::Type{<:FixedBasis}) = Base.HasLength()
Base.length(b::FixedBasis) = length(b.elts)

Base.iterate(b::FixedBasis) = iterate(b.elts)
Base.iterate(b::FixedBasis, state) = iterate(b.elts, state)
Base.IndexStyle(::Type{<:FixedBasis{T,I,V}}) where {T,I,V} = Base.IndexStyle(V)

# To break ambiguity
Base.@propagate_inbounds Base.getindex(
    b::FixedBasis{T,I},
    i::I,
) where {T,I<:Integer} = b.elts[i]
