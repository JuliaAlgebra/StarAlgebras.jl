mutable struct DiracBasis{T,I,S,M<:DiracMStructure} <: ImplicitBasis{T,I}
    object::S # any iterable
    moperation::M

    function DiracBasis{I}(itr, operation=*) where {I}
        @assert !isempty(itr)
        mstr = DiracMStructure(operation)
        return new{eltype(itr),I,typeof(itr),typeof(mstr)}(itr, mstr)
    end
end

object(db::DiracBasis) = db.object

Base.IteratorSize(::Type{<:DiracBasis{T,I,S}}) where {T,I,S} = Base.IteratorSize(S)
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

mstructure(db::DiracBasis{T}) where {T} = db.moperation

mutable struct AugmentedBasis{T,I,A<:AugmentedDirac{T},B<:AbstractBasis{T,I}} <: ImplicitBasis{A,I}
    basis::B
end

function AugmentedBasis(basis::DiracBasis{T,I}) where {T,I}
    return AugmentedBasis{T,I,AugmentedDirac{T,Int},typeof(basis)}(basis)
end

object(ab::AugmentedBasis) = object(ab.basis)

Base.IteratorSize(::Type{<:AugmentedBasis{T,A,I,B}}) where {T,A,I,B} = Base.IteratorSize(B)
function Base.length(ab::AugmentedBasis)
    @assert Base.haslength(ab.basis)
    return length(ab.basis) - 1
end

Base.iterate(ab::AugmentedBasis) = ((v, st) = iterate(object(ab)); (AugmentedDirac(v), st))
Base.iterate(ab::AugmentedBasis, st) = ((v, st) = iterate(object(ab), st); (AugmentedDirac(v), st))

Base.in(g, ab::AugmentedBasis) = false
Base.in(ad::AugmentedDirac, ab::AugmentedBasis) = ad.dirac in ab.basis

function Base.getindex(ab::AugmentedBasis{T,I,A}, x::A) where {T,I,A}
    @assert x.dirac.element in object(ab)
    return x
end

mstructure(db::AugmentedBasis) = AugmentedMStructure(mstructure(db.basis))

### translating between bases
function coeffs!(
    res::AbstractCoefficients,
    cfs::AbstractCoefficients,
    source::DiracBasis,
    target::AugmentedBasis,
)
    s = sum(values(cfs))
    if !iszero(s)
        throw("Conversion to $target not possible due to non-zero augmentation: $s")
    end
    for (k, v) in nonzero_pairs(cfs)
        isone(k) && continue
        x = source[k]
        MA.operate!(UnsafeAddMul(*), res, v, SparseCoefficients((target[AugmentedDirac(x)],), (1,)))
    end
    return MA.operate!!(canonical, res)
end
