struct ProductNotWellDefined <: Exception
    i::Any
    j::Any
    msg::Any
end

function Base.showerror(io::IO, ex::ProductNotWellDefined)
    print(
        io,
        "Product of elements $(ex.i) and $(ex.j) is not defined on the basis",
    )
    print(io, " or the multiplicative structure could not be completed")
    if isdefined(ex, :msg)
        print(io, ": $(ex.msg)")
    end
    return print(io, ".")
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

function MA.operate_to!(res, ms::MultiplicativeStructure, args::Vararg{Any,N}) where {N}
    if any(Base.Fix1(===, res), args)
        throw(ArgumentError("No alias allowed"))
    end
    MA.operate!(zero, res)
    MA.operate!(UnsafeAddMul(ms), res, args...)
    MA.operate!(canonical, res)
    return res
end

struct One end
Base.:*(::One, α) = α

function operate_with_constant!(::UnsafeAddMul, res, α, c)
    for (k, v) in nonzero_pairs(c)
        unsafe_push!(res, k, α * v)
    end
    return res
end

function operate_with_constant!(op::UnsafeAddMul, res, α, b, c, args::Vararg{Any, N}) where {N}
    for (kb, vb) in nonzero_pairs(b)
        for (kc, vc) in nonzero_pairs(c)
            operate_with_constant!(op, res, α * vb * vc, op.structure(kb, kc), args...)
        end
    end
    return res
end

_aggregate_constants(constant, non_constant) = (constant, non_constant)

function _aggregate_constants(constant, non_constant, α, args::Vararg{Any,N}) where {N}
    return _aggregate_constants(constant * α, non_constant, args...)
end

function _aggregate_constants(constant, non_constant, c::AbstractCoefficients, args::Vararg{Any,N}) where {N}
    return _aggregate_constants(constant, (non_constant..., c), args...)
end

function MA.operate!(
    op::UnsafeAddMul,
    mc::AbstractCoefficients,
    args::Vararg{Any,N},
) where {N}
    constant, non_constant = _aggregate_constants(One(), tuple(), args...)
    return operate_with_constant!(op, mc, constant, non_constant...)
end

struct DiracMStructure{Op} <: MultiplicativeStructure
    op::Op
end

function (mstr::DiracMStructure)(x::T, y::T) where {T}
    xy = mstr.op(x, y)
    return SparseCoefficients((xy,), (1,))
end
