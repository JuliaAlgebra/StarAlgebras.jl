"""
    MTable{T, I} <: MultiplicativeStructure{T}
Multiplicative table, stored explicitly as an AbstractMatrix{I}.

!!! note
    Accessing `mt[i,j]` with negative indices returns the product of
    `star`ed basis elements, i.e.
    ```julia
    mt[-i, j] == b[star(b[i])*b[j]]
    ```
"""
struct MTable{T,I,V<:AbstractVector,M<:AbstractMatrix,Ms} <:
       MultiplicativeStructure
    elts::V
    relts::Dict{T,I}
    starof::Vector{I}
    table::M
    mstr::Ms
end

function MTable(
    elts::AbstractVector,
    mstr::MultiplicativeStructure,
    dims::NTuple{2,I},
) where {I<:Integer}
    Base.require_one_based_indexing(elts)
    @assert length(elts) ≥ first(dims)
    @assert length(elts) ≥ last(dims)

    relts = Dict(b => I(idx) for (idx, b) in nonzero_pairs(elts))
    starof = [relts[star(x)] for x in elts]
    T = typeof(mstr(first(elts), first(elts)))
    table = Matrix{T}(undef, dims)
    @assert !isbitstype(T) || dims == (0, 0)

    return MTable(elts, relts, starof, table, mstr)
end

Base.@propagate_inbounds function __absindex(mt::MTable, i::Integer)
    return ifelse(i > 0, i, oftype(i, mt.starof[abs(i)]))
end

Base.size(mt::MTable, i::Vararg) = size(mt.table, i...)
Base.haskey(mt::MTable, x) = haskey(mt.relts, x)
Base.getindex(mt::MTable{T}, x::T) where {T} = mt.relts[x]
Base.@propagate_inbounds Base.getindex(mt::MTable, i::Integer) =
    mt.elts[__absindex(mt, i)]

Base.@propagate_inbounds __iscomputed(mt::MTable, i, j) =
    isassigned(mt.table, i, j) && !iszero(mt.table[i, j])

Base.@propagate_inbounds function (mt::MTable)(i::Integer, j::Integer)
    i = __absindex(mt, i)
    j = __absindex(mt, j)
    if checkbounds(Bool, mt.table, i, j)
        @inbounds begin
            if !__iscomputed(mt, i, j)
                complete!(mt, i, j)
            end
            return mt.table[i, j]
        end
    else
        g, h = mt[i], mt[j]
        return mt.mstr(g, h)
    end
end
Base.@propagate_inbounds function (mt::MTable{T})(x::T, y::T) where {T}
    i, j = mt[x], mt[y]
    return @inbounds mt(i, j)
end

Base.@propagate_inbounds function complete!(mt::MTable, i::Integer, j::Integer)
    g, h = mt[i], mt[j]
    gh_cfs = mt.mstr(g, h)
    mt.table[i, j] = gh_cfs
    return mt
end

function complete!(mt::MTable)
    Threads.@threads for j in axes(mt.table, 2)
        for i in axes(mt.table, 1)
            @inbounds complete!(mt, i, j)
        end
    end
    return mt
end

function MA.operate!(
    ms::UnsafeAddMul{<:MTable},
    res::AbstractSparseVector,
    v::AbstractVector,
    w::AbstractVector,
)
    l1 = issparse(v) ? nnz(v) : 2^8
    l2 = issparse(w) ? nnz(w) : 2^8
    l = l1 * l2
    idcs = Vector{key_type(res)}()
    vals = Vector{eltype(res)}()
    sizehint!(idcs, l)
    sizehint!(vals, l)

    for (kv, a) in nonzero_pairs(v)
        for (kw, b) in nonzero_pairs(w)
            c = ms.structure(kv, kw)
            for (k, v) in nonzero_pairs(c)
                push!(idcs, ms.structure[k])
                push!(vals, v * a * b)
            end
        end
    end
    res .+= sparsevec(idcs, vals, length(res))
    return res
end

function MA.operate!(
    ms::UnsafeAddMul{<:MTable},
    res::AbstractVector,
    v::AbstractVector,
    w::AbstractVector,
)
    for (kv, a) in nonzero_pairs(v)
        for (kw, b) in nonzero_pairs(w)
            c = ms.structure(kv, kw)
            for (k, v) in nonzero_pairs(c)
                res[ms.structure[k]] += v * a * b
            end
        end
    end
    return res
end
