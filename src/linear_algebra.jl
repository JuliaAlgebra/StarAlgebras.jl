# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

for (A, B) in (
    (Union{AlgebraElement,Term}, Union{AlgebraElement,Term}),
    (Number, Union{AlgebraElement,Term}),
    (Union{AlgebraElement,Term}, Number),
)
    @eval function MA.operate_to!(
        output::VecOrMat{<:AlgebraElement},
        ::typeof(*),
        A::AbstractMatrix{<:$A},
        B::AbstractVecOrMat{<:$B},
    )
        return _matrix_product_to!(output, A, B)
    end
end

function _matrix_product_to!(
    output::VecOrMat{P},
    A,
    B,
) where {P<:AlgebraElement}
    MA._dim_check(output, A, B)
    if Base.mightalias(output, A) || Base.mightalias(output, B)
        throw(
            ArgumentError(
                "The output of matrix multiplication must not alias an input",
            ),
        )
    end
    isempty(output) && return output
    if isempty(B)
        # No input values supply a parent for an empty contraction.
        MA.operate!(zero, output)
        return output
    end
    alg = _product_algebra(A, B)
    A = _promote_product_array(A, alg)
    B = _promote_product_array(B, alg)
    prototype = first(A)
    if first(B) isa AlgebraElement &&
       (!(prototype isa AlgebraElement) || coeffs(first(B)) isa DenseArray)
        prototype = first(B)
    elseif prototype isa Number
        prototype = first(B)
    end
    for i in eachindex(output)
        # Each output needs independent storage, including mutable zeros.
        output[i] = _sum_zero(prototype, eltype(P))
    end
    return MA.operate!(MA.add_mul, output, A, B)
end
