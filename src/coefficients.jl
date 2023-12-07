"""
    abstract type AbstractCoefficients{V,K} end
Implements `Base.keys`, `Base.values`.
"""
abstract type AbstractCoefficients{K,V} end

Base.keytype(::Type{<:AbstractCoefficients{K}}) where {K} = K
Base.valtype(::Type{<:AbstractCoefficients{K,V}}) where {K,V} = V
Base.keytype(b::AbstractCoefficients) = keytype(typeof(b))
Base.valtype(b::AbstractCoefficients) = valtype(typeof(b))

Base.iszero(sc::AbstractCoefficients) = isempty(keys(sc))

function Base.:(==)(ac1::AbstractCoefficients, ac2::AbstractCoefficients)
    for x in (ac1, ac2)
        if !__iscanonical(x)
            __canonicalize!(x)
        end
    end
    return keys(ac1) == keys(ac2) && values(ac1) == values(ac2)
end

struct DiracDelta{K,V} <: AbstractCoefficients{K,V}
    element::K
    value::V
end
DiracDelta(x) = DiracDelta(x, 1)
Base.getindex(δ::DiracDelta{K,V}, i::K) where {K,V} =
    ifelse(i == δ.element, δ.value, zero(δ.value))

__iscanonical(::DiracDelta) = true

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

function MA.operate!(::typeof(zero), s::SparseCoefficients)
    empty!(s.basis_elements)
    empty!(s.values)
    return s
end

Base.keys(sc::SparseCoefficients) = sc.basis_elements
Base.values(sc::SparseCoefficients) = sc.values

Base.zero(sc::SparseCoefficients) = SparseCoefficients(empty(keys(sc)), empty(values(sc)))

function Base.similar(s::SparseCoefficients, ::Type{T}) where {T}
    return SparseCoefficients(similar(s.basis_elements), similar(s.values, T))
end

function unsafe_append!(mc::SparseCoefficients, p::Pair{<:AbstractCoefficients,T}) where {T}
    c, val = p
    append!(mc.basis_elements, keys(c))
    for v in values(c)
        push!(mc.values, val * v)
    end
    return mc
end

function fmac!(
    res::AlgebraElement{A,T,<:SparseCoefficients},
    X::AlgebraElement{A,T,<:SparseCoefficients},
    Y::AlgebraElement{A,T,<:SparseCoefficients},
) where {A,T}
    x = coeffs(X)
    y = coeffs(Y)
    for (kx, vx) in zip(x.basis_elements, x.values)
        for (ky, vy) in zip(y.basis_elements, y.values)
            MA.operate!(MA.add_mul, res, vx * vy, basis(parent(X))[kx], basis(parent(Y))[ky])
        end
    end
    return res
end


function __canonicalize!(res::SparseCoefficients)
@error "in __canonicalize!: Not implemented yet"
    return res
end
