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
    MultiplicativeStructure{T}
Structure representing multiplication w.r.t its basis.

Implements
* `basis(ms::MultiplicativeStructure{T}) → AbstractBasis{T}`
* `Basis.getindex(ms::MultiplicativeStructure{T}, i::T, j::T) →
        Union{AbstractCoefficients, AbstractVector}`
   the product of `i` and `j` represented by coefficients in `basis(ms)`.

When the product is not representable faithfully,
   `ProductNotWellDefined` exception should be thrown.
"""
abstract type MultiplicativeStructure{T} end

function mul!(
    ms::MultiplicativeStructure,
    res::SparseCoefficients,
    v::AbstractCoefficients,
    w::AbstractCoefficients,
)
    for (kv, a) in pairs(v)
        for (kw, b) in pairs(w)
            c = ms[kv, kw] # c is an AbstractCoefficients
            unsafe_append!(res, c => a * b)
        end
    end
    __canonicalize!(res)
    return res
end

struct LazyMStructure{T,B<:AbstractBasis{T}} <: MultiplicativeStructure{T}
    basis::B
end

basis(mstr::LazyMStructure) = mstr.basis

function Base.getindex(mstr::LazyMStructure{T,<:DiracBasis{T}}, x::T, y::T) where {T}
    xy = x * y
    xy in basis(mstr) || throw(ProductNotWellDefined(x, x, "$x · $y = $xy"))
    return DiracDelta(xy)
end
