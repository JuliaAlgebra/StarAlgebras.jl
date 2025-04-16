# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

mutable struct DiracBasis{T,S,M<:DiracMStructure} <: ImplicitBasis{T,T}
    object::S # any iterable
    moperation::M

    function DiracBasis(itr, operation = *)
        @assert !isempty(itr)
        mstr = DiracMStructure(operation)
        return new{eltype(itr),typeof(itr),typeof(mstr)}(itr, mstr)
    end
end

object(db::DiracBasis) = db.object
mstructure(db::DiracBasis{T}) where {T} = db.moperation

function Base.IteratorSize(::Type{<:DiracBasis{T,S}}) where {T,S}
    return Base.IteratorSize(S)
end
function Base.length(db::DiracBasis)
    @assert Base.haslength(object(db))
    return length(object(db))
end
Base.iterate(db::DiracBasis) = iterate(object(db))
Base.iterate(db::DiracBasis, st) = iterate(object(db), st)

Base.in(g, db::DiracBasis) = g in object(db)

function Base.getindex(db::DiracBasis{T}, x::T) where {T}
    @assert x in object(db)
    return x
end
