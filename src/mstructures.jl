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

function MA.operate!(::UnsafeAddMul, res, c)
    for (k, v) in nonzero_pairs(c)
        unsafe_push!(res, k, v)
    end
    return res
end

function MA.operate!(
    op::UnsafeAddMul,
    res,
    b,
    c,
    args::Vararg{Any,N};
    cfs = nothing,
) where {N}
    for (kb, vb) in nonzero_pairs(b)
        for (kc, vc) in nonzero_pairs(c)
            for (k, v) in nonzero_pairs(op.structure(kb, kc))
                _cfs = isnothing(cfs) ? vb * vc * v : vb * vc * v * cfs
                MA.operate!(
                    op,
                    res,
                    SparseCoefficients((_key(op.structure, k),), (_cfs,)),
                    args...,
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

"""
    QuadraticForm(Q)
    QuadraticForm{star}(Q)
A simple wrapper for representing a quadratic form.

`QuadraticForm(Q)` represents an honest quadratic form, however in the context
of `*`-algebras a more canonical object is a **sesquilinear form** which can be
constructed as `QuadraticForm{star}(Q)`.

Both objects are defined by matrix coefficients with respect to their bases `b`,
i.e. their values at `i`-th and `j`-th basis elements:
    `star.(b)·Q·b = ∑ᵢⱼ Q[i,j]·star(b[i])·b[j]`.

# Necessary methods:
 * `basis(Q)` - a **finite** basis `b` of the quadratic form;
 * `Base.getindex(Q, i::T, j::T)` - the value at `i`-th and `j`-th basis elements
   where `T` is the type of indicies of `basis(Q)`;
 * `Base.eltype(Q)` - the type of `Q[i,j]`.
"""
struct QuadraticForm{T,involution}
    Q::T
    QuadraticForm(Q) = new{typeof(Q),identity}(Q)
    QuadraticForm{star}(Q) = new{typeof(Q),star}(Q)
end

Base.eltype(qf::QuadraticForm) = eltype(qf.Q)
basis(qf::QuadraticForm) = basis(qf.Q)
Base.getindex(qf::QuadraticForm, i::T, j::T) where {T} = qf.Q[i, j]

function MA.operate!(op::UnsafeAddMul, res, Q::QuadraticForm{T,ε}) where {T,ε}
    for (i, b1) in pairs(basis(Q))
        b1★ = ε(b1)
        for (j, b2) in pairs(basis(Q))
            MA.operate!(op, res, coeffs(b1★), coeffs(b2); cfs = Q[i, j])
        end
    end
    return res
end
