function _star_of(basis::AbstractBasis, len::Integer)
    return [basis[star(basis[i])] for i in 1:len]
end

"""
    MTable{I}
Multiplicative table, stored explicitly as an AbstractMatrix{I}.

Requires integer keys in `basis(mt)`.

!!! note
    Accessing `mt[i,j]` with negative indices returns the product of
    `star`ed basis elements, i.e.
    ```julia
    b = basis(mt)
    mt[-i, j] == b[star(b[i])*b[j]]
    ```
"""
struct MTable{I,M<:AbstractMatrix{I},B<:AbstractBasis{I}} <: MultiplicativeStructure{I}
    table::M
    star_of::Vector{I}
    basis::B
end

function MTable(basis::AbstractBasis{V,K}; size::Tuple{Int,Int}) where {K<:Integer}
    return MTable(zeros(K, size), _star_of(basis, max(size...)), basis)
end

basis(mt::MTable) = mt.basis
Base.size(mt::MTable) = size(mt.table)

Base.@propagate_inbounds __iscomputed(mt::MTable, i, j) = !iszero(mt.table[i, j])
function Base.@propagate_inbounds Base.getindex(mt::MTable, i::Integer, j::Integer)
    i = ifelse(i ≥ 0, i, oftype(i, mt.star_of[abs(i)]))
    j = ifelse(j ≥ 0, j, oftype(j, mt.star_of[abs(j)]))
    @boundscheck checkbounds(mt.table, i, j)

    @inbounds if !__iscomputed(mt, i, j)
        __compute_product!(mt, i, j)
    end
    @inbounds return mt.table[i, j]
end

Base.@propagate_inbounds function __compute_product!(mt::MTable, i::Integer, j::Integer)
    b = basis(mt)
    g, h = b[i], b[j]
    gh = g * h
    gh in b || throw(ProductNotWellDefined(i, j, "$g · $h = $gh"))
    mt.table[i, j] = b[gh]
    return mt
end

function complete!(table::AbstractMatrix, basis::AbstractBasis{V,K}) where {K<:Integer}
    Base.require_one_based_indexing(table)
    Threads.@threads for j in axes(table, 2)
        y = basis[j]
        for i in axes(table, 1)
            xy = basis[i] * y
            table[i, j] = basis[xy]
        end
    end
    return table
end

complete!(mt::MTable) = complete!(mt.table, basis(mt))
