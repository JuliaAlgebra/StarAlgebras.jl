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
    MultiplicativeStructure{I}
Structure representing multiplication w.r.t its basis.

Implements
* `basis(ms::MultiplicativeStructure{I}) → AbstractBasis{T,I}`
* `Basis.getindex(ms::MultiplicativeStructure{I}, i::I, j::I) →
        Union{AbstractCoefficients, AbstractVector}`
   the product of `i` and `j` represented by coefficients in `basis(ms)`.

When the product is not representable faithfully,
   `ProductNotWellDefined` exception should be thrown.
"""
abstract type MultiplicativeStructure{I} end

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

struct LazyMStructure{I,B<:AbstractBasis} <: MultiplicativeStructure{I}
    basis::B
end

LazyMStructure(basis::AbstractBasis{T,I}) where {T,I} =
    LazyMStructure{I,typeof(basis)}(basis)

basis(mstr::LazyMStructure) = mstr.basis
Base.@propagate_inbounds function Base.getindex(
    mstr::LazyMStructure{I},
    g::I,
    h::I,
) where {I}
    gh = g * h
    gh in basis(mstr) || throw(ProductNotWellDefined(i, j, "$g · $h = $gh"))
    return DiracDelta(gh)
end
