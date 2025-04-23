# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# example implementation of abstract coefficients
struct ACoeffs{T} <: SA.AbstractCoefficients{UInt32,T}
    vals::Vector{T}
end

# convienience:
ACoeffs(v::AbstractVector) = ACoeffs{valtype(v)}(v)

# AbstractCoefficients API

## Basic API
Base.keys(ac::ACoeffs) = (k for (k, v) in pairs(ac.vals) if !iszero(v))
Base.values(ac::ACoeffs) = (v for v in ac.vals if !iszero(v))
MA.operate!(::typeof(SA.canonical), ac::ACoeffs) = ac
function SA.star(b::SA.AbstractBasis, ac::ACoeffs)
    return ACoeffs([ac.vals[star(b, k)] for k in keys(ac)])
end

# for preallocation in arithmetic
function Base.similar(ac::ACoeffs, ::Type{T}) where {T}
    vals = similar(ac.vals, T)
    MA.operate!(zero, vals)
    return ACoeffs(vals)
end

# the default arithmetic implementation uses this access
Base.getindex(ac::ACoeffs, idx) = ac.vals[idx]
Base.setindex!(ac::ACoeffs, val, idx) = ac.vals[idx] = val
function SA.unsafe_push!(ac::ACoeffs, idx, val)
    ac.vals[idx] = MA.add!!(ac.vals[idx], val)
    return ac
end
