"""
    abstract type AbstractCoefficients{V,K} end
Implements `Base.keys`, `Base.values`.
"""
abstract type AbstractCoefficients{K,V} end

struct DiracDelta{K,V} <: AbstractCoefficients{K,V}
    element::K
    value::V
end
DiracDelta(x) = DiracDelta(x, 1)
Base.getindex(δ::DiracDelta{K,V}, i::K) where {K,V} =
    ifelse(i == δ.element, δ.value, zero(δ.value))

Base.keys(δ::DiracDelta) = (δ.element,)
Base.values(δ::DiracDelta) = (δ.value,)

function Base.show(io::IO, δ::DiracDelta)
    ioc = IOContext(io, :limit => true)
    print(ioc, "DiracDelta at ")
    if isone(δ.value)
        print(ioc, δ.element)
    else
        print(ioc, δ.value, '*', δ.element)
    end
end
struct SparseCoefficients{
    K,
    V,
    Vk<:AbstractVector{K},
    Vv<:AbstractVector{V}
} <: AbstractCoefficients{K,V}
    basis_elements::Vk
    values::Vv
end

Base.keys(sc::SparseCoefficients) = sc.basis_elements
Base.values(sc::SparseCoefficients) = sc.values

Base.zero(sc::SparseCoefficients) = SparseCoefficients(empty(keys(sc)), empty(values(sc)))

function unsafe_append!(mc::SparseCoefficients, p::Pair{<:AbstractCoefficients,T}) where {T}
    c, val = p
    append!(mc.basis_elements, keys(c))
    for v in values(c)
        push!(mc.values, val * v)
    end
    return mc
end

function __canonicalize!(res::SparseCoefficients)
    @error "in __canonicalize!: Not implemented yet"
    return res
end
