# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# An iterator over all natural numbers
struct NaturalNumbers end
Base.IteratorSize(::Type{NaturalNumbers}) = Base.IsInfinite()
Base.eltype(::NaturalNumbers) = Int
Base.iterate(::NaturalNumbers) = (1, 1)
Base.iterate(::NaturalNumbers, state) = (state + 1, state + 1)
Base.in(i::Int, ::NaturalNumbers) = i >= 1
Base.one(::NaturalNumbers) = 1

function Base.require_one_based_indexing(::SA.MappedBasis{T,Int,NaturalNumbers}) where {T} end
