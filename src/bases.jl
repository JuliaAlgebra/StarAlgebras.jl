"""
    AbstractBasis{T,I}
Implements a bijection between basis elements and integers.

 * `in(x, basis)` for basis elements
 * `Base.getindex(A, i::I) → T`
 * `Base.getindex(A, t::T) → I` # the bijection part

"""
abstract type AbstractBasis{T,I} end

Base.eltype(::Type{<:AbstractBasis{T}}) where {T} = T
Base.eltype(b::AbstractBasis) = eltype(typeof(b))
Base.keytype(::Type{<:AbstractBasis{T,I}}) where {T,I} = I
Base.keytype(b::AbstractBasis) = keytype(typeof(b))

"""
    ImplicitBasis{T,I}
Implicit bases are not stored in memory and can be potentially infinite.
"""
abstract type ImplicitBasis{T,I} <: AbstractBasis{T,I} end

mutable struct DiracBasis{T,I,S,St} <: ImplicitBasis{T,I}
    basis::Vector{T}
    dict::Dict{T,I}
    object::S # any iterable
    state::St
    lock::Threads.SpinLock

    function DiracBasis{I}(itr) where {I}
        elt, state = let k = iterate(itr)
            @assert !isnothing(k)
            k
        end
        T, S, St = typeof(elt), typeof(itr), typeof(state)
        return new{T,I,S,St}([elt], Dict(elt => I(1)), itr, state, Threads.SpinLock())
    end
end

function Base.length(db::DiracBasis)
    @assert Base.haslength(db)
    return length(db.object)
end

object(db::DiracBasis) = db.object
Base.eltype(db::DiracBasis) = eltype(object(db))
Base.IteratorSize(::Type{<:DiracBasis{T,I,S}}) where {T,I,S} = Base.IteratorSize(S)

Base.iterate(db::DiracBasis) = iterate(object(db))
Base.iterate(db::DiracBasis, st) = iterate(object(db), st)

Base.in(g, db::DiracBasis) = g in object(db)

function Base.getindex(db::DiracBasis{T}, x::T) where {T}
    @assert x in object(db)
    return DiracDelta(x)
end

function __cache(db::DiracBasis{T}, g::T) where {T}
    g in db || throw(KeyError(g))
    if !haskey(db.dict, g)
        lock(db.lock) do
            if !haskey(db.dict, g) # it could be computed while we waited
                k = iterate(db.object, db.state)
                while !isnothing(k)
                    h, st = k
                    if !haskey(db.dict, h)
                        push!(db.basis, h)
                        db.dict[h] = length(db.basis)
                    end
                    if h == g
                        db.state = st
                        break
                    end
                    k = iterate(db.object, db.state)
                end
            end
        end
    end
    return db.dict[g]
end

"""
    ExplicitBasis
Explicit bases are stored in an `AbstractVector` and hence immutable
(e.g. fixed in length).
"""
abstract type ExplicitBasis{T,I} <: AbstractBasis{T,I} end

function __star_of!(star_of::Vector{<:Integer}, basis::AbstractBasis{T,<:Integer}) where {T}
    for idx in eachindex(star_of)
        star_of[idx] = basis[star(basis[idx])]
    end
    return star_of
end
struct FixedBasis{T,I,A<:AbstractVector{T}} <: ExplicitBasis{T,I}
    basis::A
    rbasis::Dict{T,I}
    star_of::Vector{I}
end

function FixedBasis{I}(elts::AbstractVector) where {I}
    Base.require_one_based_indexing(elts)
    length(elts) <= typemax(I) ||
        throw("index type $I is to small for basis of length $(length(elts))")
    @assert !(eltype(elts) <: Integer)
    fb = FixedBasis(elts, Dict(b => I(idx) for (idx, b) in pairs(elts)), Vector{I}(undef, length(elts)))
    __star_of!(fb.star_of, fb)
    return fb
end

Base.size(b::FixedBasis) = size(b.basis)
Base.length(b::FixedBasis) = length(b.basis)
Base.keys(b::FixedBasis) = keys(b.basis)
Base.iterate(b::FixedBasis) = iterate(b.basis)
Base.iterate(b::FixedBasis, state) = iterate(b.basis, state)
Base.IndexStyle(::Type{<:FixedBasis{T,I,A}}) where {T,I,A} = Base.IndexStyle(A)

Base.in(g, b::FixedBasis) = haskey(b.rbasis, g)

Base.@propagate_inbounds Base.getindex(b::FixedBasis{T,I}, i::I) where {T,I} = b.basis[i]
Base.@propagate_inbounds Base.getindex(b::FixedBasis{T}, g::T) where {T} = b.rbasis[g]

# convenience only:
Base.@propagate_inbounds function Base.getindex(b::FixedBasis{T,I}, i::Integer) where {T,I<:Integer}
    i = ifelse(i > 0, i, oftype(i, b.star_of[abs(i)]))
    idx = convert(keytype(b), i)
    return b[idx]
end
# To break ambiguity
Base.@propagate_inbounds Base.getindex(b::FixedBasis{T,I}, i::I) where {T,I<:Integer} = b.basis[i]

function star(basis::FixedBasis, coeffs::SparseVector)
    nzidcs = basis.star_of[SparseArrays.nonzeroinds(coeffs)]
    nzvals = copy(SparseArrays.nonzeros(coeffs))

    v = sparsevec(nzidcs, nzvals, length(coeffs))
    dropzeros!(v)
    return v
end
