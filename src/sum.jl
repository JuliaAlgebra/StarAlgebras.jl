# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

function _sum_zero(a::AlgebraElement, ::Type{T}) where {T}
    return MA.operate!(zero, similar(a, T))
end

function _sum_zero(t::Term, ::Type{T}) where {T}
    return zero(T, parent(t))
end

# Establish a common parent for the accumulator.
function _sum_init(a, init)
    if isempty(a)
        # With no input parent, use a supplied initial value or a typed zero.
        z =
            init === nothing ?
            zero(MA.promote_operation(+, eltype(a), eltype(a))) :
            MA.copy_if_mutable(init)
        return z
    end
    alg = parent(first(a))
    T = _coeff_type(first(a))
    for x in a
        if parent(x) != alg
            alg = first(promote_bases(alg, parent(x)))
        end
        T = MA.promote_operation(+, T, _coeff_type(x))
    end
    if init !== nothing
        T = MA.promote_operation(+, T, _coeff_type(init))
        if init isa Union{AlgebraElement,Term} && parent(init) != alg
            alg = first(promote_bases(alg, parent(init)))
        end
    end
    z = _sum_zero(first(a), T)
    if parent(z) != alg
        z = first(promote_bases(z, alg))
    end
    return init === nothing ? z : z + init
end

function MA.operate(
    ::typeof(sum),
    a::AbstractArray{<:Union{AlgebraElement,Term}};
    init = nothing,
)
    z = _sum_init(a, init)
    prepare = x ->
        parent(x) == parent(z) ? x : first(promote_bases(x, parent(z)))
    return mapreduce(prepare, MA.add!!, a; init = z)
end

function Base.sum(
    a::AbstractArray{<:Union{AlgebraElement,Term}};
    dims = :,
    init = nothing,
)
    if dims === Colon()
        return MA.operate(sum, a; init)
    end
    z = _sum_init(a, init)
    # Base shares `init` across output entries, so addition must not mutate it.
    return mapreduce(identity, Base.add_sum, a; dims, init = z)
end

"""
    sum_products(a, b)

Sum the products of corresponding entries in `a` and `b`. At least one array
contains algebra elements or terms; the other may also contain numbers.
Numeric factors convert to the coefficient type of their corresponding algebra
element or term, as in ordinary scalar multiplication. Promote algebra operands
to a common basis before accumulating with MutableArithmetics. The inputs are
not modified.
"""
function sum_products end

for (A, B) in (
    (Union{AlgebraElement,Term}, Union{AlgebraElement,Term}),
    (Number, Union{AlgebraElement,Term}),
    (Union{AlgebraElement,Term}, Number),
)
    @eval function sum_products(a::AbstractArray{<:$A}, b::AbstractArray{<:$B})
        return _sum_products(a, b)
    end
end

function _sum_product_coeff_types(a, b)
    return (_coeff_type(a), _coeff_type(b))
end
function _sum_product_coeff_types(::Number, b::Union{AlgebraElement,Term})
    return (_coeff_type(b), _coeff_type(b))
end
function _sum_product_coeff_types(a::Union{AlgebraElement,Term}, ::Number)
    return (_coeff_type(a), _coeff_type(a))
end

_product_values(a) = a
_product_values(a::SparseMatrixCSC) = nonzeros(a)

function _product_prototype(a)
    values = _product_values(a)
    return isempty(values) ? nothing : first(values)
end
_product_prototype(::AbstractArray{<:Number}) = nothing

function _product_prototype(a, b)
    x, y = _product_prototype(a), _product_prototype(b)
    if y isa AlgebraElement &&
       (!(x isa AlgebraElement) || coeffs(y) isa DenseArray)
        return y
    end
    return x === nothing ? y : x
end

function _product_algebra(a, b, prototype = _product_prototype(a, b))
    alg = parent(prototype)
    for array in (a, b)
        eltype(array) <: Number && continue
        for x in _product_values(array)
            if parent(x) != alg
                alg = first(promote_bases(alg, parent(x)))
            end
        end
    end
    return alg
end

function _promote_product_array(a::AbstractArray{<:Number}, alg)
    return a
end
function _promote_product_array(a, alg)
    prepare = x -> parent(x) == alg ? x : first(promote_bases(x, alg))
    return map(prepare, a)
end

function _promote_product_array(
    a::SparseMatrixCSC{<:Union{AlgebraElement,Term}},
    alg,
)
    # Mapping a sparse matrix directly also evaluates f(zero(eltype(a))).
    # Algebra elements need a parent to construct that zero, so only map storage.
    return SparseMatrixCSC(
        size(a)...,
        copy(a.colptr),
        copy(rowvals(a)),
        _promote_product_array(nonzeros(a), alg),
    )
end

for (Wrapper, op) in
    ((LinearAlgebra.Transpose, transpose), (LinearAlgebra.Adjoint, adjoint))
    @eval begin
        _product_values(a::$Wrapper{<:Any,<:SparseMatrixCSC}) =
            $op(nonzeros(parent(a)))

        function _promote_product_array(
            a::$Wrapper{<:Union{AlgebraElement,Term},<:SparseMatrixCSC},
            alg,
        )
            # Apply the wrapper before promoting the resulting algebra elements.
            return _promote_product_array(copy(a), alg)
        end
    end
end

function _sum_products(a, b)
    MA._check_same_length(a, b)
    if isempty(a)
        return MA.fused_map_reduce(MA.add_mul, a, b)
    end
    prototype = first(a)
    if !(prototype isa AlgebraElement) &&
       first(b) isa Union{AlgebraElement,Term}
        prototype = first(b)
    end
    alg = _product_algebra(a, b)
    T = MA.promote_operation(*, _sum_product_coeff_types(first(a), first(b))...)
    for (x, y) in zip(a, b)
        T = MA.promote_operation(
            MA.add_mul,
            T,
            _sum_product_coeff_types(x, y)...,
        )
    end
    z = _sum_zero(prototype, T)
    if parent(z) != alg
        z = first(promote_bases(z, alg))
    end
    return MA.fused_map_reduce(
        MA.add_mul,
        _promote_product_array(a, alg),
        _promote_product_array(b, alg);
        init = z,
    )
end
