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

Sum the products of corresponding algebra elements or terms in `a` and `b`.
Promote all operands to a common basis before accumulating with
MutableArithmetics. The inputs are not modified.
"""
function sum_products(
    a::AbstractArray{<:Union{AlgebraElement,Term}},
    b::AbstractArray{<:Union{AlgebraElement,Term}},
)
    MA._check_same_length(a, b)
    if isempty(a)
        return MA.fused_map_reduce(MA.add_mul, a, b)
    end
    alg = parent(first(a))
    for x in Iterators.flatten((a, b))
        if parent(x) != alg
            alg = first(promote_bases(alg, parent(x)))
        end
    end
    T = MA.promote_operation(*, _coeff_type(first(a)), _coeff_type(first(b)))
    for (x, y) in zip(a, b)
        T = MA.promote_operation(MA.add_mul, T, _coeff_type(x), _coeff_type(y))
    end
    prototype = first(a) isa AlgebraElement ? first(a) : first(b)
    z = _sum_zero(prototype, T)
    if parent(z) != alg
        z = first(promote_bases(z, alg))
    end
    prepare = let alg = alg
        x -> parent(x) == alg ? x : first(promote_bases(x, alg))
    end
    return MA.fused_map_reduce(
        MA.add_mul,
        map(prepare, a),
        map(prepare, b);
        init = z,
    )
end
