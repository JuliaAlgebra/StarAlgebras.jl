abstract type AbstractMTable{I,Tw} <: MultiplicativeStructure{Tw,I} end

Base.size(mt::AbstractMTable) = size(mt.table)

_check(mt::AbstractMTable) = _check(mt.table, basis(mt), _istwisted(mt))

function _check(product_matrix, basis, twisted::Bool)
    idx = findfirst(iszero, product_matrix)
    if idx != nothing
        i, j = Tuple(idx)
        x = (twisted ? star(basis[i]) : basis[i])
        throw(
            ProductNotDefined(
                i,
                j,
                "$x · $(basis[j]) = $(_product(Val(twisted), x, basis[j]))",
            ),
        )
    end
    return true
end

_iscached(mt::AbstractMTable, i, j) = !iszero(mt.table[i, j])


## MTables

struct MTable{I,Twisted,M<:AbstractMatrix{I}} <: AbstractMTable{I,Twisted}
    table::M
end

MTable{Tw}(mt::AbstractMatrix{<:Integer}) where {Tw} = MTable{eltype(mt),Tw,typeof(mt)}(mt)

MTable(b::AbstractBasis; table_size) = MTable{false}(b; table_size=table_size)

function MTable{Tw}(basis::AbstractBasis; table_size) where {Tw}
    @assert length(table_size) == 2
    @assert 1 <= first(table_size) <= length(basis)
    @assert 1 <= last(table_size) <= length(basis)

    table = zeros(SparseArrays.indtype(basis), table_size)

    complete!(table, basis, Val(Tw))

    _check(table, basis, Tw)

    return MTable{Tw}(table)
end

for twisted in (:true, :false)
    @eval begin
        function complete!(table, basis, v::Val{$twisted})
            Threads.@threads for j in 1:size(table, 2)
                y = basis[j]
                for i in 1:size(table, 1)
                    table[i, j] = basis[_product(v, basis[i], y)]
                end
            end
            return table
        end
    end
end

basis(mt::MTable) = throw("No basis is defined for a simple $(typeof(mt))")

Base.@propagate_inbounds function Base.getindex(m::MTable, i::Integer, j::Integer)
    @boundscheck checkbounds(Bool, m, i, j) ||
                 throw(ProductNotDefined(i, j, "out of Mtable bounds"))
    @boundscheck iszero(m.table[i, j]) && throw(ProductNotDefined(i, j))
    return m.table[i, j]
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
