mutable struct FixedBasis{T,I,V<:AbstractVector{T},M<:MTable{T,I}} <:
               ExplicitBasis{T,I}
    elts::V
    table::M
end

function FixedBasis(basis::AbstractBasis; n::Integer, mt::Integer)
    @assert 0 < mt ≤ n
    elts = Iterators.take(basis, n)
    return FixedBasis(elts, mstructure(basis), (mt, mt))
end

function FixedBasis(
    elts::AbstractVector,
    mstr::MultiplicativeStructure,
    dims::NTuple{2,I},
) where {I<:Integer}
    @assert 0 < dims[1] ≤ length(elts)
    @assert 0 < dims[2] ≤ length(elts)
    @assert !(eltype(elts) <: Integer)
    return FixedBasis(elts, MTable(elts, mstr, dims))
end

mstructure(fb::FixedBasis) = fb.table

Base.in(x, b::FixedBasis) = haskey(mstructure(b), x)
Base.getindex(b::FixedBasis{T}, x::T) where {T} = mstructure(b)[x]
Base.getindex(b::FixedBasis, i::Integer) = mstructure(b)[i]

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

function coeffs(
    cfs::AbstractCoefficients,
    source::AbstractBasis,
    target::FixedBasis,
)
    source === target && return cfs
    vals = collect(values(cfs))
    idcs = collect(target[k] for k in keys(cfs))

    return sparsevec(idcs, vals, length(target))
end
