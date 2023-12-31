struct ProductNotWellDefined <: Exception
    i::Any
    j::Any
    msg::Any
end

function Base.showerror(io::IO, ex::ProductNotWellDefined)
    print(io, "Product of elements $(ex.i) and $(ex.j) is not defined on the basis")
    print(io, " or the multiplicative structure could not be completed")
    if isdefined(ex, :msg)
        print(io, ": $(ex.msg)")
    end
    print(io, ".")
end

"""
    MultiplicativeStructure
Structure representing multiplication w.r.t its basis.

Implements
* `basis(ms::MultiplicativeStructure{T}) → AbstractBasis{T}`
* `Basis.getindex(ms::MultiplicativeStructure{T}, i::T, j::T) →
        Union{AbstractCoefficients, AbstractVector}`
   the product of `i` and `j` represented by coefficients in `basis(ms)`.

When the product is not representable faithfully,
   `ProductNotWellDefined` exception should be thrown.
"""
abstract type MultiplicativeStructure end

function mul!(
    ms::MultiplicativeStructure,
    res::SparseCoefficients,
    v::AbstractCoefficients,
    w::AbstractCoefficients,
)
    MA.operate!(zero, res)
    res = fmac!(ms, res, v, w)
    __canonicalize!(res)
    return res
end

function fmac!(
    ms::MultiplicativeStructure,
    res::SparseCoefficients,
    v::AbstractCoefficients,
    w::AbstractCoefficients,
)
    for (kv, a) in pairs(v)
        for (kw, b) in pairs(w)
            c = ms(kv, kw) # ::AbstractCoefficients
            unsafe_append!(res, c => a * b)
        end
    end
    return res
end

function mul!(
    ms::MultiplicativeStructure,
    res::AbstractVector,
    X::AbstractVector,
    Y::AbstractVector,
)
    res = (res === X || res === Y) ? zero(res) : (res .= zero(eltype(res)))
    return fmac!(ms, res, X, Y)
end

struct DiracMStructure{Op} <: MultiplicativeStructure
    op::Op
end

DiracMStructure() = DiracMStructure(*)

(mstr::DiracMStructure)(x::T, y::T) where {T} = Dirac(mstr.op(x, y))
(mstr::DiracMStructure)(δx::Dirac, δy::Dirac) = mstr.op(δx, δy)

struct AugmentedMStructure{M<:DiracMStructure} <: MultiplicativeStructure
    op::M
end

function (mstr::AugmentedMStructure)(aδx::AugmentedDirac, aδy::AugmentedDirac)
    δxy = mstr.op(aδx.dirac, aδy.dirac)# :: Dirac
    c = ifelse(isone(δxy.element), zero(δxy.value), one(δxy.value))
    aδxy = AugmentedDirac(δxy)

    #(x-1)*(y-1) = 1 - x - y + xy = -1·(x-1) - 1·(y-1) + 1·(xy-1)
    return SparseCoefficients((aδx, aδy, aδxy), (-one(c), -one(c), c))
end
