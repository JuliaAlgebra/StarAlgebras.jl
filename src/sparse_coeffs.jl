struct SparseCoefficients{K,V,Vk,Vv} <: AbstractCoefficients{K,V}
    basis_elements::Vk
    values::Vv
end

function SparseCoefficients(elts::Ks, vals::Vs) where {Ks,Vs}
    return SparseCoefficients{eltype(elts),eltype(vals),Ks,Vs}(elts, vals)
end

function MA.operate!(::typeof(zero), s::SparseCoefficients)
    empty!(s.basis_elements)
    empty!(s.values)
    return s
end

Base.keys(sc::SparseCoefficients) = sc.basis_elements
Base.values(sc::SparseCoefficients) = sc.values

function Base.zero(sc::SparseCoefficients)
    return SparseCoefficients(empty(keys(sc)), empty(values(sc)))
end

function Base.similar(s::SparseCoefficients, ::Type{T}) where {T}
    return SparseCoefficients(similar(s.basis_elements), similar(s.values, T))
end

function __iscanonical(res::SparseCoefficients)
    issorted(keys(res)) || return false
    allunique(keys(res)) || return false
    return true
end

function __canonicalize!(res::SparseCoefficients)
    if !issorted(keys(res))
        p = sortperm(res.basis_elements)
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
