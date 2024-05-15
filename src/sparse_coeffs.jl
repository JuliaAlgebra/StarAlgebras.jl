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

function Base.getindex(sc::SparseCoefficients{K}, key::K) where {K}
    k = searchsortedfirst(sc.basis_elements, key; lt = comparable(K))
    if k in eachindex(sc.basis_elements)
        v = sc.values[k]
        return ifelse(sc.basis_elements[k] == key, v, zero(v))
    else
        return zero(valtype(sc))
    end
end

function Base.setindex!(sc::SparseCoefficients{K}, val, key::K) where {K}
    k = searchsortedfirst(sc.basis_elements, key; lt = comparable(K))
    if k in eachindex(sc.basis_elements) && sc.basis_elements[k] == key
        sc.values[k] += val
    else
        insert!(sc.basis_elements, k, key)
        insert!(sc.values, k, val)
    end
    return sc
end

function Base.zero(sc::SparseCoefficients)
    return SparseCoefficients(empty(keys(sc)), empty(values(sc)))
end

function Base.similar(s::SparseCoefficients, ::Type{T} = valtype(s)) where {T}
    return SparseCoefficients(similar(s.basis_elements), similar(s.values, T))
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
    T = Base._return_type(op, Tuple{valtype(X),typeof(a)})
    return similar(X, T)
end

function __prealloc(X::SparseCoefficients, Y::SparseCoefficients, op)
    # this is not even correct for op = *
    T = Base._return_type(op, Tuple{valtype(X),valtype(Y)})
    return similar(X, T)
end

comparable(::Type) = isless
function MA.operate!(::typeof(canonical), res::SparseCoefficients)
    return MA.operate!(canonical, res, comparable(key_type(res)))
end

function MA.operate!(::typeof(canonical), res::SparseCoefficients, cmp)
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
