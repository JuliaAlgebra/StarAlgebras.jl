abstract type AbstractMTable{I} <: MultiplicativeStructure{I} end

function _check(product_matrix::AbstractMatrix, basis::AbstractVector)
    idx = findfirst(iszero, product_matrix)
    if idx !== nothing
        i, j = Tuple(idx)
        x, y = basis[i], basis[j]
        throw(
            ProductNotDefined(
                i,
                j,
                "$x · $y = $(x * y)",
            ),
        )
    end
    return true
end
Base.size(mt::AbstractMTable) = size(mt.table)

function _star_of(basis::AbstractBasis, len::Integer)
    return [basis[star(basis[i])] for i in 1:len]
end

## MTables
struct MTable{I,M<:AbstractMatrix{I}} <: AbstractMTable{I}
    table::M
    star_of::Vector{I}
end

function MTable(basis::AbstractBasis; table_size)
    @assert length(table_size) == 2

    @assert 1 <= first(table_size) <= length(basis)
    @assert 1 <= last(table_size) <= length(basis)

    table = zeros(SparseArrays.indtype(basis), table_size)

    complete!(table, basis)
    _check(table, basis)

    return MTable(table, _star_of(basis, max(table_size...)))
end

function complete!(table::AbstractMatrix, basis::AbstractBasis, lck=Threads.SpinLock())
    Threads.@threads for j in axes(table, 2)
        y = basis[j]
        for i in axes(table, 1)
            xy = basis[i] * y
            lock(lck) do
                table[i, j] = xy
            end
        end
    end
    return table
end

function complete!(table::Matrix, basis::AbstractBasis, lck=Threads.SpinLock())
    Threads.@threads for j in axes(table, 2)
        y = basis[j]
        for i in axes(table, 1)
            x = basis[i]
            table[i, j] = basis[x*y]
        end
    end
    return table
end

basis(mt::MTable) = throw("No basis is defined for a simple $(typeof(mt))")
Base.@propagate_inbounds _iscached(mt::MTable, i, j) = !iszero(mt.table[i, j])
Base.@propagate_inbounds _get(cmt::MTable, i::Integer) = ifelse(i ≥ 0, i, cmt.star_of[abs(i)])

Base.@propagate_inbounds function Base.getindex(m::MTable, i::Integer, j::Integer)
    @boundscheck checkbounds(Bool, m, abs(i), abs(j)) ||
                 throw(ProductNotDefined(i, j, "out of Mtable bounds"))
    @boundscheck !_iscached(m, abs(i), abs(j)) && throw(ProductNotDefined(i, j, "product not stored"))
    @inbounds begin
        i = _get(m, i)
        j = _get(m, j)

        return m.table[i, j]
    end
end

## CachedMTables

struct CachedMTable{T,I,B<:Basis{T,I},M,Twisted} <: AbstractMTable{I,Twisted}
    basis::B
    table::M
end

CachedMTable(basis::AbstractBasis; table_size) =
    CachedMTable{false}(basis; table_size=table_size)

function CachedMTable{Tw}(basis::AbstractBasis{T,I}; table_size) where {Tw,T,I}
    return CachedMTable{Tw}(basis, zeros(I, table_size))
end

CachedMTable(basis::AbstractBasis, mt::AbstractMatrix{<:Integer}) =
    CachedMTable{false}(basis, mt)

function CachedMTable{Tw}(
    basis::AbstractBasis{T,I},
    mt::AbstractMatrix{<:Integer},
) where {Tw,T,I}
    return CachedMTable{T,I,typeof(basis),typeof(mt),Tw}(basis, mt)
end

basis(m::CachedMTable) = m.basis

Base.@propagate_inbounds function Base.getindex(cmt::CachedMTable, i::Integer, j::Integer)
    cache!(cmt, i, j)
    return cmt.table[i, j]
end

Base.@propagate_inbounds function cache!(cmt::CachedMTable, i::Integer, j::Integer)
    @boundscheck checkbounds(Bool, cmt, i, j) ||
                 throw(ProductNotDefined(i, j, "out of Mtable bounds"))
    if !_iscached(cmt, i, j)
        b = basis(cmt)
        g, h = b[i], b[j]
        gh = _product(cmt, g, h)
        gh in b || throw(ProductNotDefined(i, j, "$g · $h = $gh"))
        cmt.table[i, j] = b[gh]
    end
    return cmt
end

function cache!(
    cmt::CachedMTable,
    suppX::AbstractVector{<:Integer},
    suppY::AbstractVector{<:Integer},
)
    Threads.@threads for j in suppY
        for i in suppX
            if !_iscached(cmt, i, j)
                cache!(cmt, i, j)
            end
        end
    end
    return cmt
end

complete!(cmt::CachedMTable) = cache!(cmt, 1:size(cmt, 1), 1:size(cmt, 2))
