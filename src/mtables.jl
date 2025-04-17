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
struct MTable{T,I<:Integer,B<:AbstractBasis{T,I},Ms,M<:AbstractMatrix} <:
       MultiplicativeStructure{T,I}
    basis::B
    mstr::Ms
    starof::Vector{I}
    table::M
    lock::Base.Threads.SpinLock
end

function MTable(
    basis::AbstractBasis{T,I},
    mstr::MultiplicativeStructure,
    dims::NTuple{2,I},
) where {T,I<:Integer}
    Base.require_one_based_indexing(basis)
    @assert length(basis) ≥ first(dims)
    @assert length(basis) ≥ last(dims)

    starof = [basis[star(x)] for x in basis]
    C = typeof(mstr(first(basis), first(basis)))
    table = Matrix{C}(undef, dims)
    @assert !isbitstype(C) || dims == (0, 0)

    return MTable(basis, mstr, starof, table, Base.Threads.SpinLock())
end

Base.@propagate_inbounds function __absindex(mt::MTable, i::Integer)
    return ifelse(i > 0, i, oftype(i, mt.starof[abs(i)]))
end

Base.size(mt::MTable, i::Vararg) = size(mt.table, i...)
Base.getindex(mt::MTable, i::Integer) = mt.basis[__absindex(mt, i)]

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
