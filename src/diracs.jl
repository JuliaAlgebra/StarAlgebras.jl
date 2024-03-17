struct Dirac{K,V} <: AbstractCoefficients{K,V}
    element::K
    value::V
end
Dirac(x) = Dirac(x, 1)

Base.getindex(δ::Dirac{K,V}, i::K) where {K,V} =
    ifelse(i == δ.element, δ.value, zero(δ.value))

canonical(δ::Dirac) = δ
Base.iszero(δ::Dirac) = iszero(δ.value)

Base.keys(δ::Dirac) = (δ.element,)
Base.values(δ::Dirac) = (δ.value,)

function Base.show(io::IO, ::MIME"text/plain", δ::Dirac)
    ioc = IOContext(io, :limit => true)
    print(ioc, "DiracDelta at ")
    if isone(δ.value)
        print(ioc, δ.element)
    else
        print(ioc, δ.value, '*', δ.element)
    end
end

function Base.show(io::IO, δ::Dirac)
    ioc = IOContext(io, :limit => true)
    if isone(δ.value)
        print(ioc, δ.element)
    else
        print(ioc, δ.value, '*', δ.element)
    end
end

Base.isless(δx::Dirac, δy::Dirac) = isless(δx.element, δy.element)

# convenience:
Base.:*(δx::Dirac, δy::Dirac) = Dirac(δx.element * δy.element, δx.value * δy.value)

aug(::Dirac) = 1
