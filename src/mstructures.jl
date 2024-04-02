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

"""
    struct UnsafeAddMul{M<:Union{typeof(*),MultiplicativeStructure}}
        structure::M
    end

The value of `(op::UnsafeAddMul)(a, b, c)` is `a + structure(b, c)`
where `a` is not expected to be canonicalized before the operation `+`
and should not be expected to be canonicalized after either.
"""
struct UnsafeAddMul{M<:Union{typeof(*),MultiplicativeStructure}}
    structure::M
end

function MA.operate_to!(res, ms::MultiplicativeStructure, v, w)
    if res === v || res === w
        throw(ArgumentError("No alias allowed"))
    end
    MA.operate!(zero, res)
    res = MA.operate!(UnsafeAddMul(ms), res, v, w)
    return MA.operate!!(canonical, res)
end

function MA.operate!(
    ::UnsafeAddMul{typeof(*)},
    mc::SparseCoefficients,
    val,
    c::AbstractCoefficients,
)
    append!(mc.basis_elements, keys(c))
    vals = values(c)
    if vals isa AbstractVector
        append!(mc.values, val .* vals)
    else
        append!(mc.values, val * collect(values(c)))
    end
    return mc
end

function MA.operate!(ms::UnsafeAddMul, res, v, w)
    for (kv, a) in nonzero_pairs(v)
        for (kw, b) in nonzero_pairs(w)
            c = ms.structure(kv, kw)
            MA.operate!(UnsafeAddMul(*), res, a * b, c)
        end
    end
    return res
end

struct DiracMStructure{Op} <: MultiplicativeStructure
    op::Op
end

DiracMStructure() = DiracMStructure(*)

function (mstr::DiracMStructure)(x::T, y::T) where {T}
    xy = mstr.op(x, y)
    return SparseCoefficients((xy,), (1,))
end

struct AugmentedMStructure{M<:DiracMStructure} <: MultiplicativeStructure
    op::M
end

function (mstr::AugmentedMStructure)(aδx::Augmented, aδy::Augmented)
    δxy = first(keys(mstr.op(aδx.elt, aδy.elt)))
    if isone(δxy)
        return SparseCoefficients((aδx, aδy), (-1, -1))
    else
        aδxy = Augmented(δxy)
        #(x-1)*(y-1) = 1 - x - y + xy = -1·(x-1) - 1·(y-1) + 1·(xy-1)
        return SparseCoefficients((aδx, aδy, aδxy), (-1, -1, 1))
    end
end
