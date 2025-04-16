# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

struct SparseCoefficients{K,V,Vk,Vv} <: AbstractCoefficients{K,V}
    basis_elements::Vk
    values::Vv
end

function SparseCoefficients(elts::Ks, vals::Vs) where {Ks,Vs}
    return SparseCoefficients{eltype(elts),eltype(vals),Ks,Vs}(elts, vals)
end

Base.keys(sc::SparseCoefficients) = sc.basis_elements
Base.values(sc::SparseCoefficients) = sc.values
function Base.copy(sc::SparseCoefficients)
    return SparseCoefficients(copy(keys(sc)), copy(values(sc)))
end

function _search(keys::Tuple, key)
    # `searchsortedfirst` is not defined for `Tuple`
    return findfirst(isequal(key), keys)
end

function _search(keys, key::K) where {K}
    return searchsortedfirst(keys, key; lt = comparable(K))
end

function Base.getindex(sc::SparseCoefficients{K}, key::K) where {K}
    k = _search(sc.basis_elements, key)
    if k in eachindex(sc.basis_elements)
        v = sc.values[k]
        if sc.basis_elements[k] == key
            return v
        else
            return zero(v)
        end
    else
        return zero(value_type(sc))
    end
end

function Base.setindex!(sc::SparseCoefficients{K}, val, key::K) where {K}
    k = searchsortedfirst(sc.basis_elements, key; lt = comparable(K))
    if k in eachindex(sc.basis_elements) && sc.basis_elements[k] == key
        sc.values[k] = val
    else
        insert!(sc.basis_elements, k, key)
        insert!(sc.values, k, val)
    end
    return sc
end

################
# Broadcasting #
################

struct BroadcastStyle{K} <: Broadcast.BroadcastStyle end

Base.broadcastable(sc::SparseCoefficients) = sc
Base.BroadcastStyle(::Type{<:SparseCoefficients{K}}) where {K} = BroadcastStyle{K}()
# Disallow mixing broadcasts.
function Base.BroadcastStyle(::BroadcastStyle, ::Base.BroadcastStyle)
    return throw(
        ArgumentError(
            "Cannot broadcast `StarAlgebras.SparseCoefficients` with" *
            " another array of different type",
        ),
    )
end

# Allow broadcasting over scalars.
function Base.BroadcastStyle(style::BroadcastStyle, ::Base.Broadcast.DefaultArrayStyle{0})
    return style
end

# Used for broadcasting
Base.axes(sc::SparseCoefficients) = (sc.basis_elements,)
#Base.Broadcast.BroadcastStyle(::Type{<:SparseCoefficients{K}}) where {K} = SparseArrays.HigherOrderFns.SparseVecStyle()
#SparseArrays.HigherOrderFns.nonscalararg(::SparseCoefficients) = true

# `_get_arg` and `getindex` are inspired from `JuMP.Containers.SparseAxisArray`
_getindex(x::SparseCoefficients, index) = getindex(x, index)
_getindex(x::Any, ::Any) = x
_getindex(x::Ref, ::Any) = x[]

function _get_arg(args::Tuple, index)
    return (_getindex(first(args), index), _get_arg(Base.tail(args), index)...)
end
_get_arg(::Tuple{}, _) = ()

function Base.getindex(bc::Broadcast.Broadcasted{<:BroadcastStyle}, index)
    return bc.f(_get_arg(bc.args, index)...)
end

function Base.similar(bc::Broadcast.Broadcasted{<:BroadcastStyle}, ::Type{T}) where {T}
    return similar(_first_sparse_coeffs(bc.args...), T)
end

_first_sparse_coeffs(c::SparseCoefficients, args...) = c
_first_sparse_coeffs(_, args...) = _first_sparse_coeffs(args...)

function Base.zero(sc::SparseCoefficients)
    return SparseCoefficients(empty(keys(sc)), empty(values(sc)))
end

_similar(x::Tuple) = _similar(x, typeof(x[1]))
_similar(x::Tuple, ::Type{T}) where {T} = Vector{T}(undef, length(x))
_similar(x) = similar(x)
_similar(x, ::Type{T}) where {T} = similar(x, T)

_similar_type(::Type{<:Tuple}, ::Type{T}) where {T} = Vector{T}
_similar_type(::Type{V}, ::Type{T}) where {V,T} = similar_type(V, T)

function similar_type(::Type{SparseCoefficients{K,V,Vk,Vv}}, ::Type{T}) where {K,V,Vk,Vv,T}
    return SparseCoefficients{K,T,_similar_type(Vk, K),_similar_type(Vv, T)}
end

function Base.similar(s::SparseCoefficients, ::Type{T} = value_type(s)) where {T}
    return SparseCoefficients(collect(s.basis_elements), _similar(s.values, T))
end

function MA.mutability(
    ::Type{<:SparseCoefficients{K,V,Vk,Vv}},
    ::typeof(canonical),
    arg::Vararg{Type},
) where {K,V,Vk,Vv}
    # return MA.IsMutable()
    return MA.mutability(Vk)
end

### temporary convenience? how to handle this?
function __prealloc(X::SparseCoefficients, a::Number, op)
    T = MA.promote_operation(op, value_type(X), typeof(a))
    return similar(X, T)
end

function __prealloc(X::SparseCoefficients, Y::SparseCoefficients, op)
    # this is not even correct for op = *
    T = MA.promote_operation(op, value_type(X), value_type(Y))
    return similar(X, T)
end

comparable(::Type) = isless
function MA.operate!(::typeof(canonical), res::SparseCoefficients)
    return MA.operate!(canonical, res, comparable(key_type(res)))
end

function unsafe_push!(res::SparseCoefficients, key, value)
    push!(res.basis_elements, key)
    push!(res.values, value)
    return res
end

# `::C` is needed to force Julia specialize on the function type
# Otherwise, we get one allocation when we call `issorted`
# See https://docs.julialang.org/en/v1/manual/performance-tips/#Be-aware-of-when-Julia-avoids-specializing
function MA.operate!(::typeof(canonical), res::SparseCoefficients, cmp::C) where {C}
    sorted = issorted(res.basis_elements; lt = cmp)
    distinct = allunique(res.basis_elements)
    if sorted && distinct && !any(iszero, res.values)
        return res
    end

    if !sorted
        p = sortperm(res.basis_elements; lt = cmp)
        permute!(res.basis_elements, p)
        permute!(res.values, p)
    end

    todelete = BitSet()
    for i in firstindex(res.basis_elements):lastindex(res.basis_elements)-1
        if iszero(res.values[i])
            push!(todelete, i)
        elseif res.basis_elements[i] == res.basis_elements[i+1]
            res.values[i+1] += res.values[i]
            push!(todelete, i)
        end
    end
    if iszero(last(values(res)))
        push!(todelete, lastindex(res.basis_elements))
    end
    deleteat!(res.basis_elements, todelete)
    deleteat!(res.values, todelete)
    return res
end

# arithmetic on coefficients; performance overloads

function MA.operate!(::typeof(zero), s::SparseCoefficients)
    empty!(s.basis_elements)
    empty!(s.values)
    return s
end

function MA.operate_to!(
    res::SparseCoefficients,
    ::typeof(-),
    X::SparseCoefficients,
)
    return MA.operate_to!(res, *, X, -1)
end

function MA.operate_to!(
    res::SparseCoefficients,
    ::typeof(*),
    X::SparseCoefficients,
    a::Number,
)
    if res === X
        res.values .*= a
    else
        resize!(res.basis_elements, length(X.basis_elements))
        resize!(res.values, length(res.basis_elements))
        res.basis_elements .= X.basis_elements
        res.values .= a .* X.values
    end

    return res
end

function MA.operate_to!(
    res::SparseCoefficients,
    ::typeof(+),
    X::SparseCoefficients,
    Y::SparseCoefficients,
)
    if res === X
        append!(res.basis_elements, Y.basis_elements)
        append!(res.values, Y.values)
    elseif res === Y
        append!(res.basis_elements, X.basis_elements)
        append!(res.values, X.values)
    else
        MA.operate!(zero, res)
        append!(res.basis_elements, X.basis_elements, Y.basis_elements)
        append!(res.values, X.values, Y.values)
    end
    MA.operate!(canonical, res)
    return res
end
