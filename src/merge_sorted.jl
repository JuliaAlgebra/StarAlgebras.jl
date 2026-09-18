# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# Adapted from TypedPolynomials via MultivariatePolynomials

"""
    merge_sorted!(result, v1::AbstractVector, v2::AbstractVector; lt, combine, filter, rev=false)

In-place version of [`merge_sorted`](@ref) that writes the result into `result`.
`result` is resized to hold the merged output and may be `v1` or `v2`.
The merge writes backwards, then compacts the result after combining and filtering.
If `rev=true`, the comparison `lt` is reversed.
"""
function merge_sorted!(
    result,
    v1::AbstractVector,
    v2::AbstractVector;
    lt,
    combine,
    filter,
    rev = false,
)
    a = Iterators.Reverse(v1)
    b = Iterators.Reverse(v2)
    return _merge_sorted!(result, a, b; lt, combine, filter, rev)
end

# To make it work in case `result` is an alias for `a` or `b`, the trick
# is to merge it backward:
# https://stackoverflow.com/a/4553321
# Backwards iterators also allow lazy term products, computed once per entry.
function _merge_sorted!(
    result,
    a,
    b;
    lt,
    combine,
    filter,
    rev = false,
    by = identity,
)
    # Capture the original iterator bounds before resizing aliased storage.
    x, y = iterate(a), iterate(b)
    resize!(result, length(a) + length(b))
    i = lastindex(result)
    while x !== nothing || y !== nothing
        if y === nothing
            value = x[1]
            x = iterate(a, x[2])
        elseif x === nothing
            value = y[1]
            y = iterate(b, y[2])
        elseif by(x[1]) == by(y[1])
            value = combine(x[1], y[1])
            x, y = iterate(a, x[2]), iterate(b, y[2])
        elseif xor(lt(by(x[1]), by(y[1])), rev)
            value = y[1]
            y = iterate(b, y[2])
        else
            value = x[1]
            x = iterate(a, x[2])
        end
        if filter(value)
            result[i] = value
            i -= 1
        end
    end
    # In case some entries were equal, we allocated too much and added from the end
    # so now we need to compact things by moving them forward
    n = lastindex(result) - i
    if i >= firstindex(result)
        for j in 1:n
            result[firstindex(result)+j-1] = result[i+j]
        end
    end
    resize!(result, n)
    return result
end

"""
    first_of(a, b)

Return the first argument. Useful as a `combine` function for [`merge_sorted`](@ref)
to deduplicate without combining values.
"""
first_of(a, _) = a
MA.promote_operation(::typeof(first_of), ::Type{T}, ::Type) where {T} = T

"""
    merge_sorted(v1::AbstractVector, v2::AbstractVector; lt, combine, filter, rev=false)

Merge two sorted vectors `v1` and `v2` into a single sorted vector.
Elements are compared using `lt` (a less-than function).
When two elements are equal (neither is less than the other),
they are combined using `combine(x1, x2)`.
Elements for which `filter` returns `false` are dropped.
If `rev=true`, the comparison `lt` is reversed.
"""
function merge_sorted(
    v1::AbstractVector,
    v2::AbstractVector;
    lt,
    combine,
    filter,
    rev = false,
)
    T = MA.promote_operation(combine, eltype(v1), eltype(v2))
    result = Vector{T}(undef, length(v1) + length(v2))
    return merge_sorted!(result, v1, v2; lt, combine, filter, rev)
end

"""
    merge_sorted(a::Tuple, b::Tuple; lt, combine, filter, rev=false)

Merge two sorted tuples, removing duplicates. Returns a tuple.
If `rev=true`, the comparison `lt` is reversed.
"""
merge_sorted(::Tuple{}, ::Tuple{}; lt, combine, filter, rev = false) = ()
merge_sorted(a::Tuple, ::Tuple{}; lt, combine, filter, rev = false) = a
merge_sorted(::Tuple{}, b::Tuple; lt, combine, filter, rev = false) = b
function merge_sorted(a::Tuple, b::Tuple; lt, combine, filter, rev = false)
    x = first(a)
    y = first(b)
    if x == y
        z = combine(x, y)
        tail =
            merge_sorted(Base.tail(a), Base.tail(b); lt, combine, filter, rev)
        if filter(z)
            return (z, tail...)
        else
            return tail
        end
    elseif xor(lt(x, y), rev)
        return (x, merge_sorted(Base.tail(a), b; lt, combine, filter, rev)...)
    else
        return (y, merge_sorted(a, Base.tail(b); lt, combine, filter, rev)...)
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
    keys = merge_sorted(
        basis1.keys,
        basis2.keys;
        lt,
        combine = first_of,
        filter = _ -> true,
    )
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
