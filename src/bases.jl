abstract type AbstractBasis{T,I} <: AbstractVector{T} end

struct Basis{T,I,A<:AbstractVector{T}} <: AbstractBasis{T,I}
    basis::A
    rbasis::Dict{T,I}
end

function Basis{I}(basis::AbstractVector) where I
    length(basis) <= typemax(I) ||
        throw("index type $I is to small for basis of length $(length(basis))")
    @assert !(eltype(basis) <: Integer)
    return Basis(basis, Dict(b => I(idx) for (idx, b) in pairs(basis)))
end

Base.size(b::Basis) = size(b.basis)

Base.IndexStyle(::Type{<:Basis{T,I,A}}) where {T,I,A} = Base.IndexStyle(A)

Base.@propagate_inbounds Base.getindex(b::Basis, i::Integer) = b.basis[i]
Base.@propagate_inbounds Base.getindex(b::Basis{T}, g::T) where {T} = b.rbasis[g]

Base.in(g, b::Basis) = haskey(b.rbasis, g)

## are these really that useful?
SparseArrays.indtype(::Type{<:Basis{T,I}}) where {T,I} = I
SparseArrays.indtype(b::Basis) = SparseArrays.indtype(typeof(b))
