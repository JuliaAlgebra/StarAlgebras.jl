# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

aug(cfs::Any) = sum(values(cfs))
aug(a::AlgebraElement) = aug(coeffs(a))

function aug(ac::AbstractCoefficients)
    isempty(keys(ac)) && return zero(value_type(ac))
    return sum(c * aug(x) for (x, c) in nonzero_pairs(ac))
end

struct Augmented{K} <: AbstractCoefficients{K,Int}
    elt::K
end # corresponds to (elt - 1)

function Base.getindex(aδ::Augmented{K}, i::K) where {K}
    w = aug(aδ.elt)
    isone(aδ.elt) && return zero(w)
    i == aδ.elt && return w
    isone(i) && return -w
    return zero(w)
end

MA.operate!(::typeof(canonical), aδ::Augmented) = aδ

Base.keys(aδ::Augmented) = (k = keys(aδ.elt); (one(first(k)), first(k)))
function Base.values(aδ::Augmented)
    (e, x) = keys(aδ)
    return (aδ[e], aδ[x])
end

function star(basis::AbstractBasis, ad::Augmented)
    return Augmented(star(basis, ad.elt))
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

struct AugmentedBasis{T,I,A<:Augmented{T},B<:AbstractBasis{T,I}} <:
       ImplicitBasis{A,I}
    basis::B
end

function AugmentedBasis(basis::DiracBasis{T}) where {T}
    @assert one(object(basis)) in basis
    return AugmentedBasis{T,T,Augmented{T},typeof(basis)}(basis)
end

object(ab::AugmentedBasis) = object(ab.basis)

function Base.IteratorSize(::Type{<:AugmentedBasis{T,I,A,B}}) where {T,I,A,B}
    return Base.IteratorSize(B)
end
Base.haslength(ab::AugmentedBasis) = Base.haslength(ab.basis)

function Base.length(ab::AugmentedBasis)
    @assert Base.haslength(ab.basis)
    return length(ab.basis) - 1
end

function Base.iterate(ab::AugmentedBasis)
    (v, st) = iterate(object(ab))
    isone(v) && return iterate(ab, st)
    return Augmented(v), st
end

function Base.iterate(ab::AugmentedBasis, st)
    (v, st) = let k = iterate(object(ab), st)
        isnothing(k) && return nothing
        k
    end
    isone(v) && return iterate(ab, st)
    return Augmented(v), st
end

Base.in(g, ab::AugmentedBasis) = false
Base.in(ad::Augmented, ab::AugmentedBasis) = ad.elt in ab.basis

function Base.getindex(ab::AugmentedBasis{T,I,A}, x::A) where {T,I,A}
    @assert x.elt in object(ab)
    return x
end

mstructure(db::AugmentedBasis) = AugmentedMStructure(mstructure(db.basis))

struct AugmentedMStructure{T,I,M<:DiracMStructure{T,I}} <: MultiplicativeStructure{T,I}
    op::M
end

function (mstr::AugmentedMStructure)(aδx::Augmented, aδy::Augmented)
    δxy = first(keys(mstr.op(aδx.elt, aδy.elt)))
    if isone(δxy)
        return SparseCoefficients((aδx, aδy), (-1, -1))
    else
        aδxy = Augmented(δxy)
        #(x-1)*(y-1) = 1 - x - y + xy = -1·(x-1) - 1·(y-1) + 1·(xy-1)
        return SparseCoefficients((aδx, aδy, aδxy), (-1, -1, 1))
    end
end

function coeffs!(
    res::AbstractCoefficients,
    cfs::AbstractCoefficients,
    source::DiracBasis,
    target::AugmentedBasis,
)
    s = aug(cfs)
    if !iszero(s)
        throw(
            "Conversion to $target not possible due to non-zero augmentation: $s",
        )
    end
    for (k, v) in nonzero_pairs(cfs)
        isone(k) && continue
        x = source[k]
        MA.operate!(
            UnsafeAddMul(*),
            res,
            SparseCoefficients((target[Augmented(x)],), (v,)),
        )
    end
    MA.operate!(canonical, res)
    return res
end
