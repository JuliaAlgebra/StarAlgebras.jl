# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

import MutableArithmetics as MA

aug(cfs::Any) = sum(values(cfs))
aug(a::SA.AlgebraElement) = aug(SA.coeffs(a), SA.basis(a))

function aug(cfs, b::SA.AbstractBasis)
    SA.iszero(cfs) && return zero(SA.value_type(cfs))
    return sum(c * aug(b[x]) for (x, c) in SA.nonzero_pairs(cfs))
end

struct Augmented{K} <: SA.AbstractCoefficients{K,Int}
    elt::K
end # corresponds to (elt - 1)

function Base.getindex(aδ::Augmented{K}, i::K) where {K}
    w = aug(aδ.elt)
    isone(aδ.elt) && return zero(w)
    i == aδ.elt && return w
    isone(i) && return -w
    return zero(w)
end

MA.operate!(::typeof(SA.canonical), aδ::Augmented) = aδ

Base.keys(aδ::Augmented) = (k = keys(aδ.elt); (one(first(k)), first(k)))
function Base.values(aδ::Augmented)
    (e, x) = keys(aδ)
    return (aδ[e], aδ[x])
end

function SA.star(basis::SA.AbstractBasis, ad::Augmented)
    return Augmented(SA.star(basis, ad.elt))
end

function Base.show(io::IO, aδ::Augmented)
    ioc = IOContext(io, :limit => true)
    if iszero(aδ)
        print(ioc, "(0)")
    else
        e, x = keys(aδ)
        _, v = values(aδ)
        print(ioc, '(', -v, '·', e, '+', v, '·', x, ')')
    end
end

Base.isless(ad1::Augmented, ad2::Augmented) = isless(ad1.elt, ad2.elt)

aug(::Augmented) = 0

struct AugmentedBasis{T,I,A<:Augmented{T},B<:SA.AbstractBasis{T,I}} <:
       SA.ImplicitBasis{A,I}
    basis::B
end

function AugmentedBasis(basis::SA.DiracBasis{T}) where {T}
    @assert one(SA.object(basis)) in basis
    return AugmentedBasis{T,T,Augmented{T},typeof(basis)}(basis)
end

SA.object(ab::AugmentedBasis) = SA.object(ab.basis)

function Base.IteratorSize(::Type{<:AugmentedBasis{T,I,A,B}}) where {T,I,A,B}
    return Base.IteratorSize(B)
end
Base.haslength(ab::AugmentedBasis) = Base.haslength(ab.basis)

function Base.length(ab::AugmentedBasis)
    @assert Base.haslength(ab.basis)
    return length(ab.basis) - 1
end

function Base.iterate(ab::AugmentedBasis)
    (v, st) = iterate(SA.object(ab))
    isone(v) && return iterate(ab, st)
    return Augmented(v), st
end

function Base.iterate(ab::AugmentedBasis, st)
    (v, st) = let k = iterate(SA.object(ab), st)
        isnothing(k) && return nothing
        k
    end
    isone(v) && return iterate(ab, st)
    return Augmented(v), st
end

Base.in(g, ab::AugmentedBasis) = false
Base.in(ad::Augmented, ab::AugmentedBasis) = ad.elt in ab.basis

function Base.getindex(ab::AugmentedBasis{T,I,A}, x::A) where {T,I,A}
    @assert x.elt in SA.object(ab)
    return x
end

SA.mstructure(db::AugmentedBasis) = AugmentedMStructure(SA.mstructure(db.basis))

struct AugmentedMStructure{T,I,M<:SA.DiracMStructure{T,I}} <: SA.MultiplicativeStructure{T,I}
    op::M
end

function (mstr::AugmentedMStructure)(aδx::Augmented, aδy::Augmented)
    δxy = first(keys(mstr.op(aδx.elt, aδy.elt)))
    if isone(δxy)
        return SA.SparseCoefficients((aδx, aδy), (-1, -1))
    else
        aδxy = Augmented(δxy)
        #(x-1)*(y-1) = 1 - x - y + xy = -1·(x-1) - 1·(y-1) + 1·(xy-1)
        return SA.SparseCoefficients((aδx, aδy, aδxy), (-1, -1, 1))
    end
end

function SA.coeffs!(
    res::SA.AbstractCoefficients,
    cfs::SA.AbstractCoefficients,
    source::SA.DiracBasis,
    target::AugmentedBasis,
)
    s = aug(cfs)
    if !iszero(s)
        throw(
            "Conversion to $target not possible due to non-zero augmentation: $s",
        )
    end
    for (k, v) in SA.nonzero_pairs(cfs)
        isone(k) && continue
        x = source[k]
        MA.operate!(
            SA.UnsafeAddMul(*),
            res,
            SA.SparseCoefficients((target[Augmented(x)],), (v,)),
        )
    end
    MA.operate!(SA.canonical, res)
    return res
end
