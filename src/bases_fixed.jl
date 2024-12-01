abstract type FiniteSupportBasis{T,I} <: ExplicitBasis{T,I} end

"""
    supp(fb::FiniteSupportBasis)
Return the supporting elements of `fb` as an indexable vector
"""
function supp end

Base.IteratorSize(::Type{<:FiniteSupportBasis}) = Base.HasLength()
Base.IteratorEltype(::Type{<:FiniteSupportBasis{T}}) where {T} = T
Base.length(b::FiniteSupportBasis) = length(supp(b))

Base.iterate(b::FiniteSupportBasis) = iterate(supp(b))
Base.iterate(b::FiniteSupportBasis, state) = iterate(supp(b), state)
# function Base.IndexStyle(::Type{<:FiniteSupportBasis{T,I,V}}) where {T,I,V}
#     return Base.IndexStyle(V)
# end

Base.@propagate_inbounds function Base.getindex(
    b::FiniteSupportBasis,
    i::Integer,
)
    return supp(b)[i]
end

mutable struct FixedBasis{T,I,V<:AbstractVector{T},M<:MTable{T,I}} <:
               FiniteSupportBasis{T,I}
    supporting_elts::V
    table::M
end

function FixedBasis(basis::AbstractBasis; n::Integer, mt::Integer = UInt32(0))
    @assert 0 ≤ mt ≤ n
    elts = Iterators.take(basis, n)
    return FixedBasis(collect(elts), mstructure(basis), (mt, mt))
end

function FixedBasis(
    elts::AbstractVector,
    mstr::MultiplicativeStructure,
    dims::NTuple{2,I} = (UInt32(0), UInt32(0)),
) where {I<:Integer}
    @assert 0 ≤ dims[1] ≤ length(elts)
    @assert 0 ≤ dims[2] ≤ length(elts)
    @assert !(eltype(elts) <: Integer)
    return FixedBasis(elts, MTable(elts, mstr, dims))
end

supp(fb::FiniteSupportBasis) = fb.supporting_elts
mstructure(fb::FixedBasis) = fb.table
Base.in(x, b::FixedBasis) = haskey(mstructure(b), x)
Base.getindex(b::FixedBasis{T}, x::T) where {T} = mstructure(b)[x]

struct SubBasis{T,I,V<:AbstractVector{I},B<:AbstractBasis{T,I}} <:
       FiniteSupportBasis{T,I}
    supporting_idcs::V
    parent_basis::B
end

supp(sb::SubBasis) = sb.supporting_idcs
Base.parent(sub::SubBasis) = sub.parent_basis

Base.in(x, b::SubBasis) = x in supp(b)
function Base.getindex(b::SubBasis{T,I}, x::T) where {T,I<:Integer}
    return convert(I, parent(b)[supp(b)[x]])
end

function Base.getindex(b::SubBasis{T,T}, x::T) where {T}
    parent(b)[x] in supp(b) && return x
end

struct SubMStructure{SB<:SubBasis,MS} <: MultiplicativeStructure
    basis::SB
    mstructure::MS
end

function mstructure(b::SubBasis)
    pb = b.parent_basis
    return SubMStructure(b, mstructure(pb))
end

function (mstr::SubMStructure)(x::T, y::T) where {T}
    b = mstr.basis
    xy = mstr.mstructure(b[x], b[y])
    return xy
end

# this is used for key-lookup in mstructures.jl
# MA.operate!(op::UnsafeAddMul, …)
_key(mstr::SubMStructure, k) = findfirst(==(k), supp(mstr.basis))
