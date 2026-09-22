# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

for (A, B) in (
    (Union{AlgebraElement,Term}, Union{AlgebraElement,Term}),
    (Number, Union{AlgebraElement,Term}),
    (Union{AlgebraElement,Term}, Number),
)
    @eval function MA.operate_to!(
        output::Union{
            VecOrMat{<:AlgebraElement},
            SparseMatrixCSC{<:AlgebraElement},
        },
        ::typeof(*),
        A::AbstractMatrix{<:$A},
        B::AbstractVecOrMat{<:$B},
        α::Number = true,
    )
        return _matrix_product_to!(output, A, B, α)
    end
end

_matrix_product_mightalias(a, b) = Base.mightalias(a, b)

function _matrix_product_mightalias(a::SparseMatrixCSC, b::SparseMatrixCSC)
    # Distinct empty buffers can share dataids. Compare the arrays themselves
    # as well: sharing an empty, resizable buffer still aliases its owner.
    return any(
        x === y || Base.mightalias(x, y) for
        x in (a.colptr, rowvals(a), nonzeros(a)),
        y in (b.colptr, rowvals(b), nonzeros(b))
    )
end

function _matrix_product_mightalias(
    a,
    b::Union{LinearAlgebra.Transpose,LinearAlgebra.Adjoint},
)
    return _matrix_product_mightalias(a, parent(b))
end

function _check_matrix_product(output, A, B)
    MA._dim_check(output, A, B)
    if _matrix_product_mightalias(output, A) ||
       _matrix_product_mightalias(output, B)
        throw(
            ArgumentError(
                "The output of matrix multiplication must not alias an input",
            ),
        )
    end
    return
end

function _matrix_product_zeros!(output::VecOrMat{P}, prototype) where {P}
    for i in eachindex(output)
        # Each output needs independent storage, including mutable zeros.
        output[i] = _sum_zero(prototype, eltype(P))
    end
    return output
end

function _matrix_product_zeros!(output::SparseMatrixCSC, prototype)
    return MA.operate!(zero, output)
end

function _matrix_product_to!(
    output::Union{VecOrMat{P},SparseMatrixCSC{P}},
    A,
    B,
    α,
) where {P<:AlgebraElement}
    _check_matrix_product(output, A, B)
    isempty(output) && return output
    prototype = _product_prototype(A, B)
    if prototype === nothing
        # Empty storage may leave no input values supplying a parent.
        MA.operate!(zero, output)
        return output
    end
    alg = _product_algebra(A, B, prototype)
    prototype = first(promote_bases(prototype, alg))
    _matrix_product_zeros!(output, prototype)
    iszero(α) && return output
    if output isa SparseMatrixCSC &&
       A isa
       Union{SparseMatrixCSC,MA._TransposeOrAdjoint{<:Any,<:SparseMatrixCSC}} &&
       B isa
       Union{SparseMatrixCSC,MA._TransposeOrAdjoint{<:Any,<:SparseMatrixCSC}}
        init =
            (a, b) ->
                MA.operate!(MA.add_mul, _sum_zero(prototype, eltype(P)), a, b)
        MA._spmatmul!(init, output, A, B)
    else
        MA.operate!(MA.add_mul, output, A, B)
    end
    if !isone(α)
        for p in _product_values(output)
            MA.operate_to!(p, *, p, α)
        end
    end
    return output
end

function MA._mul!(
    output::VecOrMat{<:AlgebraElement},
    A::AbstractMatrix,
    B::AbstractVecOrMat,
    α::Number,
    β::Number,
)
    _check_matrix_product(output, A, B)
    isempty(output) && return output
    if iszero(β)
        # The destination may be uninitialized when its old value is unused.
        return MA.operate_to!(output, *, A, B, α)
    end
    if iszero(α) || isempty(B)
        if !isone(β)
            for i in eachindex(output)
                output[i] = output[i] * β
            end
        end
        return output
    end
    # Keep the old destination while the existing kernel computes A * B * α.
    # A reusable per-entry accumulator in MA could avoid this temporary array;
    # reusing the current matrix kernel is sufficient for now.
    product = MA.operate_to!(similar(output), *, A, B, α)
    alg = _product_algebra(product, output)
    for i in eachindex(output, product)
        p = first(promote_bases(product[i], alg))
        c = output[i]
        output[i] =
            isone(β) ? MA.operate!(+, p, c) : MA.operate!(MA.add_mul, p, c, β)
    end
    return output
end
