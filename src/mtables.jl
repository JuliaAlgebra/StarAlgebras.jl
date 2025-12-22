# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

"""
    MTable{T, I} <: MultiplicativeStructure{T,I}
Multiplicative table, stored explicitly as an AbstractMatrix{I}.

!!! note
    Accessing `mt(i,j)` with negative indices returns the product of
    `star`ed basis elements, i.e.
    ```julia
    mt(-i, j) == b[star(b[i])*b[j]]
    ```
"""
struct MTable{T,I<:Integer,Ms<:MultiplicativeStructure{T,I},M<:AbstractMatrix} <:
       MultiplicativeStructure{T,I}
    mstr::Ms
    table::M
    lock::Base.Threads.SpinLock
end

function MTable(
    mstr::MultiplicativeStructure{T,I},
    dims::NTuple{2,I},
) where {T,I<:Integer}
    bas = basis(mstr)
    Base.require_one_based_indexing(bas)

    if Base.haslength(bas)
        @assert length(bas) ≥ first(dims)
        @assert length(bas) ≥ last(dims)
    end

    C = typeof(mstr(first(bas), first(bas)))
    table = Matrix{C}(undef, dims)
    # this is to avoid situation with allocated garbage in table
    # we want table to consist of #undefs as "sentiel values"
    @assert !isbitstype(C) || dims == (0, 0)

    return MTable(mstr, table, Base.Threads.SpinLock())
end

function MTable(
    basis::AbstractBasis{T,I},
    dims::NTuple{2,I},
) where {T,I<:Integer}
    return MTable(DiracMStructure(basis, *), dims)
end

Base.@propagate_inbounds function __absindex(mt::MTable, i::Integer)
    return ifelse(i > 0, i, oftype(i, star(mt, abs(i))))
end

Base.size(mt::MTable, i::Vararg) = size(mt.table, i...)
Base.getindex(mt::MTable, i::Integer) = mt.mstr[__absindex(mt, i)]
basis(mstr::MTable) = basis(mstr.mstr)

function __iscomputed(mt::MTable, i, j)
    return isassigned(mt.table, i, j)
end

function (mt::MTable{T,I})(x, y, ::Type{I}) where {T,I}
    return map_keys(Base.Fix1(getindex, mt), mt(x, y, T))
end

(mt::MTable{T,I})(x::Integer, y::Integer) where {T,I} = mt(x, y, I)

function (mt::MTable{T})(x::T, y::T, ::Type{U}) where {T,U}
    i = __absindex(mt, mt[x])
    j = __absindex(mt, mt[y])
    return mt(i, j, U)
end

function (mt::MTable{T})(i::Integer, j::Integer, ::Type{T}) where {T}
    if !checkbounds(Bool, mt.table, i, j)
        x, y = mt[i], mt[j]
        return mt.mstr(x, y)
    end
    if !__iscomputed(mt, i, j)
        lock(mt.lock) do
            return thread_unsafe_complete!(mt, i, j)
        end
    end

    # TODO remove for v1.12 https://github.com/JuliaLang/julia/pull/54707
    res = mt.table[i, j] # load
    while res != mt.table[i, j] # compare
        res = mt.table[i, j] # and load again
        yield()
    end

    return mt.table[i, j] # and load again
end

function thread_unsafe_complete!(mt::MTable, i::Integer, j::Integer)
    g, h = mt[i], mt[j]
    gh_cfs = mt.mstr(g, h)
    # TODO use atomic for v1.12 https://github.com/JuliaLang/julia/pull/54707
    #      but then we'll have to use `AtomicMemory`
    mt.table[i, j] = gh_cfs
    return mt
end

function complete!(mt::MTable)
    Threads.@threads for j in axes(mt.table, 2)
        for i in axes(mt.table, 1)
            @inbounds thread_unsafe_complete!(mt, i, j)
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
