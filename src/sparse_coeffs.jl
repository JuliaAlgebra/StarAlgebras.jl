# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

struct SparseCoefficients{K,V,Vk,Vv,L} <: AbstractCoefficients{K,V}
    basis_elements::Vk
    values::Vv
    isless::L
end

function SparseCoefficients(elts::Ks, vals::Vs, isless = isless) where {Ks,Vs}
    return SparseCoefficients{eltype(elts),eltype(vals),Ks,Vs,typeof(isless)}(
        elts,
        vals,
        isless,
    )
end

Base.keys(sc::SparseCoefficients) = sc.basis_elements
Base.values(sc::SparseCoefficients) = sc.values
function Base.copy(sc::SparseCoefficients)
    return SparseCoefficients(copy(keys(sc)), copy(values(sc)), sc.isless)
end

function MA.mutability(
    ::Type{<:SparseCoefficients{K,V,Vector{K},Vector{V}}},
) where {K,V}
    return MA.IsMutable()
end

function MA.promote_operation(
    ::typeof(+),
    ::Type{X},
    ::Type{Y},
) where {X<:SparseCoefficients,Y<:SparseCoefficients}
    return similar_type(
        X,
        MA.promote_operation(+, value_type(X), value_type(Y)),
    )
end

function MA.mutable_copy(sc::SparseCoefficients)
    result = copy(sc)
    map!(MA.copy_if_mutable, keys(result), keys(result))
    map!(MA.copy_if_mutable, values(result), values(result))
    return result
end

function _search(keys::Tuple, key; lt)
    # `searchsortedfirst` is not defined for `Tuple`
    return findfirst(isequal(key), keys)
end

function _search(keys, key::K; lt) where {K}
    return searchsortedfirst(keys, key; lt)
end

function Base.getindex(sc::SparseCoefficients{K}, key::K) where {K}
    k = _search(sc.basis_elements, key; lt = sc.isless)
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
    k = searchsortedfirst(sc.basis_elements, key; lt = sc.isless)
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
function Base.BroadcastStyle(::Type{<:SparseCoefficients{K}}) where {K}
    return BroadcastStyle{K}()
end
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
function Base.BroadcastStyle(
    style::BroadcastStyle,
    ::Base.Broadcast.DefaultArrayStyle{0},
)
    return style
end

# Used for broadcasting
Base.axes(sc::SparseCoefficients) = (sc.basis_elements,)
#Base.Broadcast.BroadcastStyle(::Type{<:SparseCoefficients{K}}) where {K} = SparseArrays.HigherOrderFns.SparseVecStyle()
#SparseArrays.HigherOrderFns.nonscalararg(::SparseCoefficients) = true

# `_get_arg` and `_getindex` are inspired from `JuMP.Containers.SparseAxisArray`.
# We do not implement `Base.getindex` for `Broadcasted` and use our own `copy`
# and `copyto!` instead as `getindex(::Broadcasted{<:BroadcastStyle}, key)`
# would invalidate `getindex(::Broadcasted, ::Any)` calls compiled by `Base`
# and `LinearAlgebra`.
_getindex(x::SparseCoefficients, key) = getindex(x, key)
_getindex(x::Any, ::Any) = x
_getindex(x::Ref, ::Any) = x[]
function _getindex(bc::Broadcast.Broadcasted, key)
    return bc.f(_get_arg(bc.args, key)...)
end

function _get_arg(args::Tuple, key)
    return (_getindex(first(args), key), _get_arg(Base.tail(args), key)...)
end
_get_arg(::Tuple{}, _) = ()

function Base.copy(bc::Broadcast.Broadcasted{<:BroadcastStyle})
    keys = only(axes(bc))
    return SparseCoefficients(
        collect(keys),
        [_getindex(bc, key) for key in keys],
        _first_sparse_coeffs(bc.args...).isless,
    )
end

function Base.copyto!(
    dest::SparseCoefficients,
    bc::Broadcast.Broadcasted{<:BroadcastStyle},
)
    for key in only(axes(bc))
        dest[key] = _getindex(bc, key)
    end
    return dest
end

_first_sparse_coeffs(c::SparseCoefficients, args...) = c
function _first_sparse_coeffs(bc::Broadcast.Broadcasted, args...)
    return _first_sparse_coeffs(bc.args..., args...)
end
_first_sparse_coeffs(_, args...) = _first_sparse_coeffs(args...)

function Base.zero(sc::SparseCoefficients)
    return SparseCoefficients(empty(keys(sc)), empty(values(sc)), sc.isless)
end

_similar(x::Tuple) = _similar(x, typeof(x[1]))
_similar(x::Tuple, ::Type{T}) where {T} = Vector{T}(undef, length(x))
_similar(x) = similar(x)
_similar(x, ::Type{T}) where {T} = similar(x, T)

_similar_type(::Type{<:Tuple}, ::Type{T}) where {T} = Vector{T}
_similar_type(::Type{V}, ::Type{T}) where {V,T} = similar_type(V, T)

function similar_type(
    ::Type{SparseCoefficients{K,V,Vk,Vv,L}},
    ::Type{T},
) where {K,V,Vk,Vv,T,L}
    return SparseCoefficients{K,T,_similar_type(Vk, K),_similar_type(Vv, T),L}
end

function Base.similar(
    s::SparseCoefficients,
    ::Type{T} = value_type(s),
) where {T}
    return SparseCoefficients(
        collect(s.basis_elements),
        _similar(s.values, T),
        s.isless,
    )
end

function Base.convert(
    ::Type{SparseCoefficients{K,V,Vh,Vv,L}},
    s::SparseCoefficients{K,V,Vh,Vv,L},
) where {K,V,Vh,Vv,L}
    return s
end

_collect_if_needed(::Type, v) = v
_collect_if_needed(::Type{<:AbstractVector}, t::Tuple) = collect(t)

function Base.convert(
    ::Type{SparseCoefficients{K,V,Vh,Vv,L}},
    s::SparseCoefficients,
) where {K,V,Vh,Vv,L}
    return SparseCoefficients{K,V,Vh,Vv,L}(
        _collect_if_needed(Vh, s.basis_elements),
        _collect_if_needed(Vv, s.values),
        s.isless,
    )
end

function map_keys(f::Function, s::SparseCoefficients)
    return SparseCoefficients(map(f, s.basis_elements), s.values, s.isless)
end

"""
    map_keys!(f, s::SparseCoefficients)

Replace each key of `s` by `f(key)`, keeping the key storage and coefficients.
The map must preserve key order and distinctness; no sorting or merging is
performed. To preserve shared keys, `f` must not mutate its argument.
"""
function map_keys!(f::F, s::SparseCoefficients) where {F}
    map!(f, keys(s), keys(s))
    return s
end

function MA.mutability(
    ::Type{<:SparseCoefficients{K,V,Vk,Vv}},
    ::typeof(canonical),
    arg::Vararg{Type},
) where {K,V,Vk,Vv}
    # return MA.IsMutable()
    return MA.mutability(Vk)
end

function __prealloc(X::SparseCoefficients, Y::SparseCoefficients, op)
    # this is not even correct for op = *
    T = MA.promote_operation(op, value_type(X), value_type(Y))
    return similar(X, T)
end

function unsafe_push!(res::SparseCoefficients, key, value)
    push!(res.basis_elements, key)
    push!(res.values, value)
    return res
end

function MA.operate!(
    ::typeof(remove_leading_term),
    c::SparseCoefficients{K,V,Vector{K},Vector{V}},
) where {K,V}
    if !isempty(keys(c))
        pop!(c.basis_elements)
        pop!(c.values)
    end
    return c
end

# `{...,L}` is needed to force Julia specialize on the function type
# Otherwise, we get one allocation when we call `issorted`
# See https://docs.julialang.org/en/v1/manual/performance-tips/#Be-aware-of-when-Julia-avoids-specializing
function MA.operate!(
    ::typeof(canonical),
    res::SparseCoefficients{K,V,Vk,Vv,L},
) where {K,V,Vk,Vv,L}
    sorted = issorted(res.basis_elements; lt = res.isless)
    distinct = allunique(res.basis_elements)
    if sorted && distinct && !any(iszero, res.values)
        return res
    end

    if !sorted
        p = sortperm(res.basis_elements; lt = res.isless)
        permute!(res.basis_elements, p)
        permute!(res.values, p)
    end

    todelete = BitSet()
    for i in firstindex(res.basis_elements):(lastindex(res.basis_elements)-1)
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

# An indexed view of the parallel arrays for the shared merge kernel.
struct _CoefficientPairs{K,V,C<:AbstractCoefficients{K,V}} <:
       AbstractVector{Pair{K,V}}
    coefficients::C
end
Base.size(p::_CoefficientPairs) = (length(keys(p.coefficients)),)
Base.firstindex(p::_CoefficientPairs) = firstindex(keys(p.coefficients))
Base.lastindex(p::_CoefficientPairs) = lastindex(keys(p.coefficients))
function Base.getindex(p::_CoefficientPairs, i::Int)
    return keys(p.coefficients)[i] => values(p.coefficients)[i]
end
function Base.setindex!(p::_CoefficientPairs, value, i::Int)
    keys(p.coefficients)[i], values(p.coefficients)[i] = value
    return p
end
function Base.resize!(p::_CoefficientPairs, n::Int)
    resize!(keys(p.coefficients), n)
    resize!(values(p.coefficients), n)
    return p
end

function Base.deleteat!(p::_CoefficientPairs, indices)
    deleteat!(keys(p.coefficients), indices)
    deleteat!(values(p.coefficients), indices)
    return p
end

function map_coefficients_to!(
    res::SparseCoefficients,
    f::F,
    X::SparseCoefficients;
    nonzero = false,
) where {F}
    output = _CoefficientPairs(res)
    resize!(output, length(keys(X)))
    map!(p -> first(p) => f(last(p)), output, _CoefficientPairs(X))
    if !nonzero
        filter!(p -> !iszero(last(p)), output)
    end
    if !_strictly_sorted(res)
        MA.operate!(canonical, res)
    end
    return res
end

# The result owns its values so that modifying them in place, e.g., with
# `MA.operate!(*, value, 2)`, does not modify the inputs.
struct _CopyValue
    skip::Bool
end
function (c::_CopyValue)(p::Pair)
    return c.skip ? p : first(p) => MA.copy_if_mutable(last(p))
end

# Entries of `Y` that are not combined are mapped by `only_b` after
# `transform`, pass `identity` if `transform` already returns owned values.
function _merge_coefficients!(
    res,
    X,
    Y,
    transform::F = identity,
    only_b::G = _CopyValue(res === Y),
) where {F,G}
    _merge_sorted!(
        _CoefficientPairs(res),
        _CoefficientPairs(X),
        _CoefficientPairs(Y);
        lt = res.isless,
        by = first,
        map_b = transform,
        combine = (a, b) -> first(a) => last(a) + last(b),
        filter = p -> !iszero(last(p)),
        only_a = _CopyValue(res === X),
        only_b,
    )
    return res
end

# Read-only coefficients returned by a multiplicative structure can contain
# unsorted or repeated keys. Preserve the canonicalizing addition for those.
function _strictly_sorted(c::SparseCoefficients)
    k = keys(c)
    return all(i -> c.isless(k[i], k[i+1]), firstindex(k):(lastindex(k)-1))
end

function MA.operate!(::typeof(+), X::SparseCoefficients, Y::SparseCoefficients)
    return MA.operate_to!(X, +, X, Y)
end

function MA.operate!(::typeof(zero), s::SparseCoefficients)
    empty!(s.basis_elements)
    empty!(s.values)
    return s
end

function MA.operate_to!(
    res::SparseCoefficients,
    ::typeof(+),
    X::SparseCoefficients,
    Y::SparseCoefficients,
)
    if res.isless == X.isless == Y.isless &&
       _strictly_sorted(X) &&
       _strictly_sorted(Y)
        return _merge_coefficients!(res, X, Y)
    end
    copied(c) = Iterators.map(MA.copy_if_mutable, c.values)
    if res === X
        append!(res.basis_elements, Y.basis_elements)
        append!(res.values, copied(Y))
    elseif res === Y
        append!(res.basis_elements, X.basis_elements)
        append!(res.values, copied(X))
    else
        MA.operate!(zero, res)
        append!(res.basis_elements, X.basis_elements, Y.basis_elements)
        append!(res.values, copied(X), copied(Y))
    end
    MA.operate!(canonical, res)
    return res
end
