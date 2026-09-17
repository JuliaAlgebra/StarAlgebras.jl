# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

"""
    abstract type AbstractCoefficients{K,V} end
Everything that implements a fixed set of methods can be used as
`SparseCoefficients` without subtyping it.

# "Read-only" coefficients
E.g. returned by calls to a `MultiplicativeStructure` need to implement

* `Base.keys`
* `Base.values`
* `StarAlgebras.canonical`
* `StarAlgebras.star(b::AbstractBasis, ac)`

For general types either all necessary arithmetic operations need to be
implemented, or fallbacks using the framework of `MutableArithmetics` are
provided based on random indexing. Additionally one needs to provide:

* `Base.similar(ac, T::Type)` with the same semantics as the one for vectors
* `Base.getindex(ac, idx)`
* `Base.setindex!(ac, val, idx)`
* `MutableArithmetics.operate!(ms::UnsafeAddMul, ac, v::C, w::C) where C<:SA.AbstractCoefficients`

"""
abstract type AbstractCoefficients{K,V} end

key_type(::Type{<:AbstractCoefficients{K}}) where {K} = K
value_type(::Type{<:AbstractCoefficients{K,V}}) where {K,V} = V
key_type(b::AbstractCoefficients) = key_type(typeof(b))
value_type(b::AbstractCoefficients) = value_type(typeof(b))

key_type(b) = keytype(b)
# `keytype(::Type{SparseVector{V,K}})` is not defined so it falls
# back to `keytype{::Type{<:AbstractArray})` which returns `Int`.
key_type(::Type{SparseArrays.SparseVector{V,K}}) where {V,K} = K
key_type(v::SparseArrays.SparseVector) = key_type(typeof(v))
key_type(::Tuple) = Int

value_type(coeffs) = valtype(coeffs)
# `valtype` is not defined for `Tuple` so we need to define
# our own function `value_type` as defining `valtype` for
# tuples would be type piracy
value_type(::Type{NTuple{N,T}}) where {N,T} = T
value_type(::NTuple{N,T}) where {N,T} = T

Base.iszero(ac::AbstractCoefficients) = isempty(keys(ac))

Base.similar(ac::AbstractCoefficients) = similar(ac, value_type(ac))

similar_type(::Type{<:Vector}, ::Type{T}) where {T} = Vector{T}
function similar_type(
    ::Type{<:SparseArrays.SparseVector{C,I}},
    ::Type{T},
) where {C,I,T}
    return SparseArrays.SparseVector{T,I}
end

"""
    canonical(ac::AbstractCoefficients)
Compute the canonical form of `ac` (e.g. grouping coefficients together, etc).

If `ac` can be brought to canonical form in-place one has to implement
* `MA.mutability(::Type{typeof(ac)}, canonical, ::Vararg{Type}) = MA.IsMutable()`
* `MA.operate!(canonical, ac)` that performs this canonicalization.

otherwise `canonical(ac)` needs to be implemented.
"""
function canonical end
function MA.promote_operation(::typeof(canonical), ::Type{C}) where {C}
    return C
end

# example implementation for vectors
MA.operate!(::typeof(canonical), sv::SparseVector) = dropzeros!(sv)
MA.operate!(::typeof(canonical), v::Vector) = v

function Base.:(==)(ac1::AbstractCoefficients, ac2::AbstractCoefficients)
    MA.operate!(canonical, ac1)
    MA.operate!(canonical, ac2)
    length(keys(ac1)) == length(keys(ac2)) || return false
    all(x -> ==(x...), zip(keys(ac1), keys(ac2))) || return false
    all(x -> ==(x...), zip(values(ac1), values(ac2))) || return false
    return true
end

function Base.hash(ac::AbstractCoefficients, h::UInt)
    MA.operate!(canonical, ac)
    return foldl((h, i) -> hash(i, h), nonzero_pairs(ac); init = h)
end

"""
    nonzero_pairs(ac::AbstractCoefficients)
Return an iterator over pairs `(k=>v)` of keys and values stored in `ac`.

The iterator contains all pairs with `v` potentially non-zero.
"""
@inline function nonzero_pairs(ac)
    return (k => v for (k, v) in zip(keys(ac), values(ac)))
end

@inline nonzero_pairs(v::AbstractVector) =
    (p for p in pairs(v) if !iszero(last(p)))
@inline function nonzero_pairs(v::AbstractSparseVector)
    return zip(SparseArrays.nonzeroinds(v), SparseArrays.nonzeros(v))
end

function LinearAlgebra.norm(sc::AbstractCoefficients, p::Real)
    isempty(values(sc)) && return (0^p)^(1 / p)
    return sum(abs(v)^p for v in values(sc))^(1 / p)
end

function LinearAlgebra.dot(ac::AbstractCoefficients, bc::AbstractCoefficients)
    if isempty(values(ac)) || isempty(values(bc))
        return zero(MA.promote_sum_mul(value_type(ac), value_type(bc)))
    else
        return sum(c * star(bc[i]) for (i, c) in nonzero_pairs(ac))
    end
end

function LinearAlgebra.dot(w::AbstractVector, ac::AbstractCoefficients)
    @assert key_type(ac) <: Integer
    if isempty(values(ac))
        return zero(MA.promote_sum_mul(eltype(w), value_type(ac)))
    else
        return sum(w[i] * v for (i, v) in nonzero_pairs(ac))
    end
end

function LinearAlgebra.dot(ac::AbstractCoefficients, w::AbstractVector)
    @assert key_type(ac) <: Integer
    if isempty(values(ac))
        return zero(MA.promote_sum_mul(eltype(w), value_type(ac)))
    else
        return sum(v * star(w[i]) for (i, v) in nonzero_pairs(ac))
    end
end

Base.zero(X::AbstractCoefficients) = MA.operate!(zero, similar(X))
function Base.:-(X::AbstractCoefficients)
    return MA.operate_to!(
        similar(X, MA.promote_operation(-, value_type(X))),
        -,
        X,
    )
end
function Base.:*(a::Union{T,Number}, X::AbstractCoefficients{K,T}) where {K,T}
    res = similar(X, MA.promote_operation(*, T, T))
    return MA.operate_to!(res, *, convert(T, a), X)
end
for op in (:*, :/, ://, :div)
    @eval function Base.$op(
        X::AbstractCoefficients{K,T},
        a::Union{T,Number},
    ) where {K,T}
        R =
            $op === (//) ? Base.promote_op($op, T, T) :
            MA.promote_operation($op, T, T)
        res = similar(X, R)
        return MA.operate_to!(res, $op, X, convert(T, a))
    end
end
function Base.:+(X::AbstractCoefficients, Y::AbstractCoefficients)
    return MA.operate_to!(__prealloc(X, Y, +), +, X, Y)
end
function Base.:-(X::AbstractCoefficients, Y::AbstractCoefficients)
    return MA.operate_to!(__prealloc(X, Y, -), -, X, Y)
end

## fallbacks; require Base.setindex! implemented for AbstractCoefficients

function MA.operate!(::typeof(zero), X::AbstractCoefficients)
    for (idx, x) in nonzero_pairs(X)
        X[idx] = zero(x)
    end
    return X
end

# unary -
function MA.operate_to!(
    res::AbstractCoefficients,
    ::typeof(-),
    X::AbstractCoefficients,
)
    if res !== X
        MA.operate!(zero, res)
    end
    for (idx, x) in nonzero_pairs(X)
        res[idx] = -x
    end
    return res
end

# scalar
function MA.operate_to!(
    res::AbstractCoefficients,
    ::typeof(*),
    a::Union{T,Number},
    X::AbstractCoefficients{K,T},
) where {K,T}
    return _map_coefficients!(res, Base.Fix1(*, convert(T, a)), X)
end

function MA.operate_to!(
    res::AbstractCoefficients,
    op::Union{typeof(*),typeof(/),typeof(//),typeof(div)},
    X::AbstractCoefficients{K,T},
    a::Union{T,Number},
) where {K,T}
    return _map_coefficients!(res, Base.Fix2(op, convert(T, a)), X)
end

function _map_coefficients!(res, f, X)
    if res !== X
        MA.operate!(zero, res)
    end
    for (idx, x) in nonzero_pairs(X)
        res[idx] = f(x)
    end
    MA.operate!(canonical, res)
    return res
end

# binary ops
function MA.operate_to!(
    res::AbstractCoefficients,
    ::typeof(+),
    X::AbstractCoefficients,
    Y::AbstractCoefficients,
)
    if res === X
        for (idx, y) in nonzero_pairs(Y)
            res[idx] += y
        end
    elseif res === Y
        for (idx, x) in nonzero_pairs(X)
            res[idx] += x
        end
    else
        MA.operate!(zero, res)
        for (idx, x) in nonzero_pairs(X)
            res[idx] += x
        end
        for (idx, y) in nonzero_pairs(Y)
            res[idx] += y
        end
    end
    return res
end

function MA.operate_to!(
    res::AbstractCoefficients,
    ::typeof(-),
    X::AbstractCoefficients,
    Y::AbstractCoefficients,
)
    if res === X
        for (idx, y) in nonzero_pairs(Y)
            X[idx] -= y
        end
    else
        if res !== Y
            MA.operate!(zero, res)
        end
        MA.operate_to!(res, -, Y)
        MA.operate_to!(res, +, res, X)
    end
    return res
end

__lmul(a, k) = a * k
__rmul(a, k) = k * a

function MA.operate_to!(
    res::AbstractCoefficients,
    mul::Union{typeof(__lmul),typeof(__rmul)},
    a,
    X::AbstractCoefficients,
)
    if res === X
        throw("No aliasing is allowed in shift")
    else
        MA.operate!(zero, res)
        for (k, v) in nonzero_pairs(X)
            res[mul(a, k)] += v
        end
    end
    return res
end
