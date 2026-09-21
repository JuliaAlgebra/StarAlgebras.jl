# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

for (A, B) in (
    (Union{AlgebraElement,Term}, Union{AlgebraElement,Term}),
    (Number, Union{AlgebraElement,Term}),
    (Union{AlgebraElement,Term}, Number),
)
    @eval function MA.operate_to!(
        output::Vector{<:AlgebraElement},
        ::typeof(*),
        A::AbstractMatrix{<:$A},
        b::AbstractVector{<:$B},
    )
        return _matrix_vector_to!(output, A, b)
    end
end

function _matrix_vector_to!(output::Vector{P}, A, b) where {P<:AlgebraElement}
    MA._dim_check(output, A, b)
    if Base.mightalias(output, A) || Base.mightalias(output, b)
        throw(
            ArgumentError(
                "The output of matrix multiplication must not alias an input",
            ),
        )
    end
    isempty(output) && return output
    if isempty(b)
        # No input values supply a parent for an empty contraction.
        MA.operate!(zero, output)
        return output
    end
    alg = _product_algebra(A, b)
    A = _promote_product_array(A, alg)
    b = _promote_product_array(b, alg)
    prototype = first(A)
    if first(b) isa AlgebraElement &&
       (!(prototype isa AlgebraElement) || coeffs(first(b)) isa DenseArray)
        prototype = first(b)
    elseif prototype isa Number
        prototype = first(b)
    end
    for i in eachindex(output)
        # Each output needs independent storage, including mutable zeros.
        output[i] = _sum_zero(prototype, eltype(P))
    end
    return MA.operate!(MA.add_mul, output, A, b)
end
