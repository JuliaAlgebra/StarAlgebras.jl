# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

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
    MultiplicativeStructure{T,I}
Structure representing multiplication w.r.t its basis.

Implements
* `basis(ms::MultiplicativeStructure{T,I}) → AbstractBasis{T}`
* `(ms::MultiplicativeStructure{T,I})(i::V, j::V, ::Type{U}) where {V<:Union{T,I},U<:Union{T,I}} →
        Union{AbstractCoefficients, AbstractVector}`
   the product of `i` and `j` represented by coefficients in `basis(ms)`
   for which the keys are of type `U<:Union{T,I}`.

The following shortcuts are provided:
* `(ms::MultiplicativeStructure{T})(i::T, j::T)` → `ms(i, j, T)`
* `(ms::MultiplicativeStructure{T,I})(i::I, j::I)` → `ms(i, j, I)`
* `(ms::MultiplicativeStructure)[x]` → `basis(ms)[x]`

When the product is not representable faithfully,
   `ProductNotWellDefined` exception should be thrown.
"""
abstract type MultiplicativeStructure{T,I} end

function (mstr::MultiplicativeStructure{T})(x::T, y::T) where {T}
    return mstr(x, y, T)
end

function (mstr::MultiplicativeStructure{T,I})(x::Integer, y::Integer) where {T,I<:Integer}
    return mstr(convert(I, x), convert(I, y), I)
end

# To break ambiguity for `MappedBasis{T,T,S,identity,identity}` for which `T` and `I` are the same
function (mstr::MultiplicativeStructure{T,T})(x::T, y::T) where {T}
    return mstr(x, y, T)
end

Base.getindex(mstr::MultiplicativeStructure, x) = basis(mstr)[x]

basis(mstr::MultiplicativeStructure) = mstr.basis

struct UnsafeAddMul{M<:Union{typeof(*),MultiplicativeStructure}}
    structure::M
end

"""
    operate_to!(res, ms::MultiplicativeStructure, A, B, α)
Compute `α·A·B` storing the result in `res`. Return `res`.

`res` is assumed to behave like `AbstractCoefficients` and not aliased with any
other arguments.
`A` and `B` are assumed to behave like `AbstractCoefficients`, while `α` will
be treated as a scalar.

Canonicalization of the result happens only once at the end of the operation.
"""
function MA.operate_to!(res, ms::MultiplicativeStructure, A, B, α = true)
    if any(Base.Fix1(===, res), (A, B))
        throw(
            ArgumentError(
                "Aliasing arguments in multiplication is not supported",
            ),
        )
    end
    MA.operate!(zero, res)
    res = MA.operate!(UnsafeAddMul(ms), res, A, B, α)
    MA.operate!(canonical, res)
    return res
end

struct UnsafeAdd end

function MA.operate!(::UnsafeAdd, res, b)
    for (k, v) in nonzero_pairs(b)
        unsafe_push!(res, k, v)
    end
    return res
end

function MA.operate!(::UnsafeAddMul, res, A, α)
    for (k, v) in nonzero_pairs(A)
        unsafe_push!(res, k, v * α)
    end
end

function MA.operate!(op::UnsafeAddMul, res, A, B, α)
    for (kA, vA) in nonzero_pairs(A)
        for (kB, vB) in nonzero_pairs(B)
            MA.operate!(
                op,
                res,
                op.structure(kA, kB),
                MA.@rewrite(α * vA * vB),
            )
        end
    end
    return res
end

struct DiracMStructure{T,I,B<:AbstractBasis{T,I},Op} <: MultiplicativeStructure{T,I}
    basis::B
    op::Op
end

function MA.promote_operation(::typeof(basis), ::Type{<:DiracMStructure{T,I,B}}) where {T,I,B}
    return B
end

function (mstr::DiracMStructure{T})(x::T, y::T, ::Type{T}) where {T}
    xy = mstr.op(x, y)
    return SparseCoefficients((xy,), (1,))
end

function (mt::DiracMStructure{T,I})(x, y, ::Type{I}) where {T,I}
    return map_keys(Base.Fix1(getindex, mt), mt(x, y, T))
end

function (mstr::DiracMStructure{T,I})(x::I, y::I, ::Type{U}) where {T,I,U}
    return mstr(mstr[x], mstr[y], U)
end

# To break ambiguity for `MappedBasis{T,T,S,identity,identity}` for which `T` and `I` are the same
function (mstr::DiracMStructure{T,T})(x::T, y::T, ::Type{T}) where {T}
    xy = mstr.op(x, y)
    return SparseCoefficients((xy,), (1,))
end
