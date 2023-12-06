"""
    abstract type AbstractCoefficients{V,K} end
Implements `Base.keys`, `Base.values`.
"""
abstract type AbstractCoefficients{V,K} end

struct DiracDelta{V,K} <: AbstractCoefficients{V,K}
    value::V
    element::K
end
DiracDelta(x) = DiracDelta(1, x)
Base.getindex(δ::DiracDelta{V,K}, i::K) where {V,K} =
    ifelse(i == δ.element, δ.value, zero(δ.value))

Base.keys(δ::DiracDelta) = (δ.element,)
Base.values(δ::DiracDelta) = (δ.value,)

struct SparseCoefficients{V,K,Vv<:AbstractVector{V},Vk<:AbstractVector{K}} <: AbstractCoefficients{V,K}
    values::Vv
    basis_elements::Vk
end

Base.keys(sc::SparseCoefficients) = sc.basis_elements
Base.values(sc::SparseCoefficients) = sc.values

function mul!(
    ms::MultiplicativeStructure,
    res::SparseCoefficients,
    v::AbstractCoefficients,
    w::AbstractCoefficients,
)
    for (a, kv) in pairs(v)
        for (b, kw) in pairs(w)
            c = ms[kv, kw]
            unsafe_append!(res, c => a * b)
        end
    end
    __canonicalize!(res)
    return res
end

function unsafe_append!(mc::SparseCoefficients, p::Pair{<:SparseCoefficients,T}) where {T}
    c, val = p
    append!(mc.basis_elements, keys(c))
    for v in values(c)
        push!(mc.values, val * v)
    end
    return mc
end
