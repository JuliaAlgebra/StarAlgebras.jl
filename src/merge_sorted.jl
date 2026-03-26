# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# Adapted from TypedPolynomials via MultivariatePolynomials

"""
    merge_sorted!(result, v1::AbstractVector, v2::AbstractVector, lt, combine, filter)

In-place version of [`merge_sorted`](@ref) that writes the result into `result`.
`result` must be large enough to hold the merged output (at most `length(v1) + length(v2)`).
It is resized to the actual output length before returning.
"""
function merge_sorted!(
    result,
    v1::AbstractVector,
    v2::AbstractVector;
    lt,
    combine,
    filter,
)
    i = firstindex(result)
    i1 = firstindex(v1)
    i2 = firstindex(v2)
    while i1 <= lastindex(v1) && i2 <= lastindex(v2)
        x1 = v1[i1]
        x2 = v2[i2]
        if lt(x1, x2)
            if filter(x1)
                result[i] = x1
                i += 1
            end
            i1 += 1
        elseif lt(x2, x1)
            if filter(x2)
                result[i] = x2
                i += 1
            end
            i2 += 1
        else
            c = combine(x1, x2)
            if filter(c)
                result[i] = c
                i += 1
            end
            i1 += 1
            i2 += 1
        end
    end
    for j = i1:lastindex(v1)
        if filter(v1[j])
            result[i] = v1[j]
            i += 1
        end
    end
    for j = i2:lastindex(v2)
        if filter(v2[j])
            result[i] = v2[j]
            i += 1
        end
    end
    resize!(result, i - 1)
    return result
end

# TODO add this promote_operation in MA, we cannot here as it would be type piracy
_promote_op(::typeof(first), ::Type{T}, ::Type) where {T} = T
_promote_op(op, ::Type{T}, ::Type{U}) where {T,U} = MA.promote_operation(op, T, U)

"""
    merge_sorted(v1::AbstractVector, v2::AbstractVector; lt, combine, filter)

Merge two sorted vectors `v1` and `v2` into a single sorted vector.
Elements are compared using `lt` (a less-than function).
When two elements are equal (neither is less than the other),
they are combined using `combine(x1, x2)`.
Elements for which `filter` returns `false` are dropped.
"""
function merge_sorted(
    v1::AbstractVector,
    v2::AbstractVector;
    lt,
    combine,
    filter,
)
    T = _promote_op(combine, eltype(v1), eltype(v2))
    result = Vector{T}(undef, length(v1) + length(v2))
    return merge_sorted!(result, v1, v2; lt, combine, filter)
end

"""
    merge_sorted(a::Tuple, b::Tuple; lt = isless)

Merge two sorted tuples, removing duplicates. Returns a tuple.
"""
merge_sorted(::Tuple{}, ::Tuple{}; lt, combine, filter) = ()
merge_sorted(a::Tuple, ::Tuple{}; lt, combine, filter) = a
merge_sorted(::Tuple{}, b::Tuple; lt, combine, filter) = b
function merge_sorted(a::Tuple, b::Tuple; lt, combine, filter)
    x = first(a)
    y = first(b)
    if x == y
        z = combine(x, y)
        tail = merge_sorted(Base.tail(a), Base.tail(b); lt, combine, filter)
        if filter(z)
            return (z, tail...)
        else
            return tail
        end
    elseif lt(x, y)
        return (x, merge_sorted(Base.tail(a), b; lt, combine, filter)...)
    else
        return (y, merge_sorted(a, Base.tail(b); lt, combine, filter)...)
    end
end

"""
    multi_findsorted(x, y; lt = isless)

Given two sorted collections `x` and `y`, return a vector `I` of length `length(x)`
where `I[i]` is the index of `x[i]` in `y`, or `0` if `x[i]` does not occur in `y`.
Both `x` and `y` must be sorted according to `lt`.
"""
function multi_findsorted(x, y; lt = isless)
    I = zeros(Int, length(x))
    j = 1
    for i in eachindex(x)
        while j ≤ length(y) && lt(y[j], x[i])
            j += 1
        end
        if j ≤ length(y) && x[i] == y[j]
            I[i] = j
        end
    end
    return I
end

"""
    merge_bases_with_maps(basis1::SubBasis, basis2::SubBasis)

Merge two sorted `SubBasis` with the same parent basis into a single `SubBasis`.
Returns `(merged_basis, I1, I2)` where `I1` and `I2` are index vectors
such that `I1[i]` is the index of the `i`th element of `merged_basis` in `basis1`
(or `0` if absent), and similarly for `I2`.
"""
function merge_bases_with_maps(basis1::SB, basis2::SB) where {SB<:SubBasis}
    @assert basis1.parent_basis == basis2.parent_basis
    @assert basis1.is_sorted
    @assert basis2.is_sorted
    lt = comparable(parent(basis1))
    keys = merge_sorted(basis1.keys, basis2.keys; lt, combine = first, filter = _ -> true)
    I1 = multi_findsorted(keys, basis1.keys; lt)
    I2 = multi_findsorted(keys, basis2.keys; lt)
    return SubBasis(basis1.parent_basis, keys), I1, I2
end

"""
    merge_bases(basis1::SubBasis, basis2::SubBasis)

Merge two sorted `SubBasis` with the same parent basis into a single `SubBasis`.
See also [`merge_bases_with_maps`](@ref).
"""
function merge_bases(basis1::SB, basis2::SB) where {SB<:SubBasis}
    return first(merge_bases_with_maps(basis1, basis2))
end
