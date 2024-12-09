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

struct UnsafeAddMul{M<:Union{typeof(*),MultiplicativeStructure}}
    structure::M
end

"""
    operate_to!(res, ms::MultiplicativeStructure, A, B[, α = true])
Compute `α·A·B` storing the result in `res`. Return `res`.

`res` is assumed to behave like `AbstractCoefficients` and not aliased with any
other arguments.
`A` and `B` are assumed to behave like `AbstractCoefficients`, while `α` will
be treated as a scalar.

Canonicalization of the result happens only once at the end of the operation.
"""
function MA.operate_to!(res, ms::MultiplicativeStructure, A, B, α = true)
    if res === A || res === B
        throw(
            ArgumentError(
                "Aliasing arguments in multiplication is not supported",
            ),
        )
    end
    MA.operate!(zero, res)
    res = MA.operate_to!(res, UnsafeAddMul(ms), A, B, α)
    MA.operate!(canonical, res)
    return res
end

struct UnsafeAdd end

function MA.operate_to!(res, ::UnsafeAdd, b)
    for (k, v) in nonzero_pairs(b)
        unsafe_push!(res, k, v)
    end
    return res
end

function MA.operate_to!(res, op::UnsafeAddMul, A, B, α = true)
    for (kA, vA) in nonzero_pairs(A)
        for (kB, vB) in nonzero_pairs(B)
            for (k, v) in nonzero_pairs(op.structure(kA, kB))
                cfs = MA.@rewrite α * vA * vB * v
                MA.operate_to!(
                    res,
                    UnsafeAdd(),
                    SparseCoefficients((_key(op.structure, k),), (cfs,)),
                )
            end
        end
    end
    return res
end

struct DiracMStructure{Op} <: MultiplicativeStructure
    op::Op
end

function (mstr::DiracMStructure)(x::T, y::T) where {T}
    xy = mstr.op(x, y)
    return SparseCoefficients((xy,), (1,))
end
