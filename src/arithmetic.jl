# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

_coeff_type(::Type{A}) where {A<:AlgebraElement} = eltype(A)
_coeff_type(::Type{<:Term{T}}) where {T} = T
_coeff_type(a::Type) = a
_coeff_type(a) = _coeff_type(typeof(a))

function _product_coefficient_type(
    ::Type{<:Union{AlgebraElement{T,A},Term{T,A}}},
    ::Type{Y},
) where {T,A,Y}
    M = MA.promote_operation(mstructure, A)
    C = MA.promote_operation(*, T, _coeff_type(Y))
    return MA.promote_operation(*, value_type(M), C)
end

function algebra_promote_operation(op, args::Vararg{Type,N}) where {N}
    T = if op === (*) && args[1] <: AlgebraElement && args[2] <: AlgebraElement
        _product_coefficient_type(args...)
    else
        MA.promote_operation(op, _coeff_type.(args)...)
    end
    if args[2] <: AlgebraElement &&
       MA.promote_operation(coeffs, args[2]) <: DenseArray # what a hack :)
        return similar_type(args[2], T)
    end
    return similar_type(args[1], T)
end

# module structure:

function Base.:-(X::AlgebraElement)
    return MA.operate_to!(similar(X, MA.promote_operation(-, eltype(X))), -, X)
end
function MA.promote_operation(
    ::typeof(*),
    ::Type{S},
    ::Type{A},
) where {T,A<:AlgebraElement{T},S<:Union{T,Number}}
    return similar_type(A, MA.promote_operation(*, S, T))
end
function MA.promote_operation(
    op::Union{typeof(*),typeof(/),typeof(//),typeof(div)},
    ::Type{A},
    ::Type{S},
) where {T,A<:AlgebraElement{T},S<:Union{T,Number}}
    if op === (//)
        return similar_type(A, Base.promote_op(op, T, S))
    end
    return algebra_promote_operation(op, A, S)
end
function _scalar_lmul(a, X::AlgebraElement)
    P = MA.promote_operation(*, typeof(a), typeof(X))
    return MA.operate_to!(MA.zero_for(P, (X,)), *, a, X)
end
function _scalar_right(op, X::AlgebraElement, a)
    R = eltype(MA.promote_operation(op, typeof(X), typeof(a)))
    return MA.operate_to!(similar(X, R), op, X, a)
end
Base.:*(a::Union{T,Number}, X::AlgebraElement{T}) where {T} = _scalar_lmul(a, X)
for op in (:*, :/, ://, :div)
    @eval function Base.$op(X::AlgebraElement{T}, a::Union{T,Number}) where {T}
        return _scalar_right($op, X, a)
    end
end
function Base.:*(
    a::T,
    X::AlgebraElement{C,A},
) where {T,C,O,A<:AbstractStarAlgebra{O,T}}
    return MA.operate_to!(similar(X), __lmul, a, X)
end
function Base.:*(
    X::AlgebraElement{C,A},
    a::T,
) where {C,T,O,A<:AbstractStarAlgebra{O,T}}
    return MA.operate_to!(similar(X), __rmul, a, X)
end

for op in [:+, :-, :*]
    @eval begin
        function MA.promote_operation(
            ::typeof($op),
            ::Type{X},
            ::Type{Y},
        ) where {X<:AlgebraElement,Y<:AlgebraElement}
            return algebra_promote_operation($op, X, Y)
        end
        function Base.$op(X::AlgebraElement, Y::AlgebraElement)
            _X, _Y = promote_bases(X, Y)
            P = MA.promote_operation($op, typeof(_X), typeof(_Y))
            return MA.operate_to!(MA.zero_for(P, (_X,), (_Y,)), $op, _X, _Y)
        end
    end
end

Base.:^(a::AlgebraElement, p::Integer) = Base.power_by_squaring(a, p)

# The destination's parent cannot change during an in-place operation.
function _assert_same_basis(
    op,
    A::AlgebraElement,
    B::Union{AlgebraElement,Term},
)
    parent(A) == parent(B) && return
    return throw(
        ArgumentError(
            "cannot `$op` in place with this destination basis. Prepare the " *
            "destination in a common basis with `StarAlgebras.promote_bases`, " *
            "or use the `*`, `+`, `-` operators which promote automatically.",
        ),
    )
end

_promote_operand(op, output, x) = x

# Promote operands without changing the prepared destination's parent.
function _promote_operand(op, output, x::Union{AlgebraElement,Term})
    parent(x) == parent(output) && return x
    x = first(promote_bases(x, parent(output)))
    _assert_same_basis(op, output, x)
    return x
end

# mutable API

function map_coefficients_to!(
    res::AlgebraElement,
    f::F,
    X::AlgebraElement;
    nonzero = false,
) where {F}
    X = _promote_operand(map_coefficients_to!, res, X)
    map_coefficients_to!(coeffs(res), f, coeffs(X); nonzero)
    return res
end

function map_coefficients_to!(
    res::AlgebraElement,
    f::F,
    t::Term;
    nonzero = false,
) where {F}
    t = _promote_operand(map_coefficients_to!, res, t)
    if iszero(t)
        return MA.operate!(zero, res)
    end
    single =
        SparseCoefficients((MA.copy_if_mutable(t.index),), (coefficient(t),))
    map_coefficients_to!(coeffs(res), f, single; nonzero)
    return res
end

function MA.operate!(::typeof(zero), a::AlgebraElement)
    MA.operate!(zero, coeffs(a))
    return a
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(*),
    a::Union{T,Number},
    X::AlgebraElement{T},
) where {T}
    X = _promote_operand(*, res, X)
    MA.operate_to!(coeffs(res), *, a, coeffs(X))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    op::Union{typeof(*),typeof(/),typeof(//),typeof(div)},
    X::AlgebraElement{T},
    a::Union{T,Number},
) where {T}
    X = _promote_operand(op, res, X)
    MA.operate_to!(coeffs(res), op, coeffs(X), a)
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    mul::Union{typeof(__lmul),typeof(__rmul)},
    a,
    X::AlgebraElement,
)
    X = _promote_operand(mul, res, X)
    MA.operate_to!(coeffs(res), mul, a, coeffs(X))
    return res
end

function MA.operate_to!(res::AlgebraElement, ::typeof(-), X::AlgebraElement)
    X = _promote_operand(-, res, X)
    MA.operate_to!(coeffs(res), -, coeffs(X))
    return res
end

function MA.operate!(::typeof(+), X::AlgebraElement, Y::AlgebraElement)
    return MA.operate_to!(X, +, X, Y)
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(+),
    X::AlgebraElement,
    Y::AlgebraElement,
)
    X = _promote_operand(+, res, X)
    Y = _promote_operand(+, res, Y)
    MA.operate_to!(coeffs(res), +, coeffs(X), coeffs(Y))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(-),
    X::AlgebraElement,
    Y::AlgebraElement,
)
    X = _promote_operand(-, res, X)
    Y = _promote_operand(-, res, Y)
    MA.operate_to!(coeffs(res), -, coeffs(X), coeffs(Y))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(*),
    A::AlgebraElement,
    B::AlgebraElement,
)
    A = _promote_operand(*, res, A)
    B = _promote_operand(*, res, B)
    mstr = mstructure(res)
    MA.operate_to!(coeffs(res), mstr, coeffs(A), coeffs(B), true)
    return res
end

function MA.operate!(
    op::MA.AddSubMul,
    f::AlgebraElement,
    a::AlgebraElement,
    b::AlgebraElement,
)
    a = _promote_operand(op, f, a)
    b = _promote_operand(op, f, b)
    c, ca, cb = coeffs(f), coeffs(a), coeffs(b)
    if c === ca || c === cb
        # Preserve both factors before accumulation changes their storage.
        original = MA.mutable_copy(c)
        ca = c === ca ? original : ca
        cb = c === cb ? original : cb
    end
    α = op === MA.add_mul ? true : -1
    MA.operate!(UnsafeAddMul(mstructure(f)), c, ca, cb, α)
    MA.operate!(canonical, c)
    return f
end

function MA.operate_to!(
    output::AlgebraElement,
    op::MA.AddSubMul,
    f::AlgebraElement,
    a::AlgebraElement,
    b::AlgebraElement,
)
    a, b = _prepare_fused_output!(output, op, f, a, b)
    return MA.operate!(op, output, a, b)
end

function _prepare_fused_output!(output, op, f, a, b)
    f = _promote_operand(op, output, f)
    a = _promote_operand(op, output, a)
    b = _promote_operand(op, output, b)
    if coeffs(output) !== coeffs(f)
        if (a isa AlgebraElement && coeffs(output) === coeffs(a)) ||
           (b isa AlgebraElement && coeffs(output) === coeffs(b))
            throw(
                ArgumentError(
                    "Aliasing the product in operate_to! is not supported; use operate!",
                ),
            )
        end
        MA.operate!(zero, output)
        for (k, v) in nonzero_pairs(coeffs(f))
            unsafe_push!(
                coeffs(output),
                MA.copy_if_mutable(k),
                MA.copy_if_mutable(v),
            )
        end
        MA.operate!(canonical, coeffs(output))
    end
    return a, b
end

# Scalars use the same coefficient arithmetic as ordinary multiplication.
for (L, R, left) in (
    (:(Union{T,Number}), :(Union{AlgebraElement{T},Term{T}}), true),
    (:(Union{AlgebraElement{T},Term{T}}), :(Union{T,Number}), false),
)
    a, g = left ? (:a, :b) : (:b, :a)
    fix = left ? Base.Fix1 : Base.Fix2
    @eval begin
        function MA.promote_operation(
            op::MA.AddSubMul,
            ::Type{F},
            ::Type{A},
            ::Type{B},
        ) where {T,F<:AlgebraElement,A<:$L,B<:$R}
            return similar_type(
                F,
                MA.promote_operation(
                    op,
                    eltype(F),
                    _coeff_type(A),
                    _coeff_type(B),
                ),
            )
        end

        function MA.operate!(
            op::MA.AddSubMul,
            f::AlgebraElement,
            a::$L,
            b::$R,
        ) where {T}
            return _scalar_add_mul!(op, f, $g, $fix(*, $a))
        end

        function MA.operate_to!(
            output::AlgebraElement,
            op::MA.AddSubMul,
            f::AlgebraElement,
            a::$L,
            b::$R,
        ) where {T}
            a, b = _prepare_fused_output!(output, op, f, a, b)
            return MA.operate!(op, output, a, b)
        end
    end
end

function _scalar_add_mul!(
    op,
    f::AlgebraElement,
    g::AlgebraElement,
    scale::F,
) where {F}
    g = _promote_operand(op, f, g)
    _add_scaled_coefficients!(MA.add_sub_op(op), coeffs(f), coeffs(g), scale)
    return f
end

function _scalar_add_mul!(op, f::AlgebraElement, t::Term, scale::F) where {F}
    value = MA.add_sub_op(op)(scale(coefficient(t)))
    return MA.operate!(+, f, Term(parent(t), t.index, value))
end

function _add_scaled_coefficients!(add, f, g, scale::F) where {F}
    for (k, v) in nonzero_pairs(g)
        # Replace coefficients: they may be shared with the scalar or input.
        f[k] = add(f[k], scale(v))
    end
    MA.operate!(canonical, f)
    return f
end

function _add_scaled_coefficients!(
    add,
    f::SparseCoefficients,
    g::SparseCoefficients,
    scale::F,
) where {F}
    product = p -> MA.copy_if_mutable(first(p)) => add(scale(last(p)))
    if f.isless == g.isless && _strictly_sorted(f) && _strictly_sorted(g)
        return _merge_coefficients!(f, f, g, product)
    end
    # Capture the original range so this also works when f === g.
    for i in eachindex(keys(g))
        k, v = product(keys(g)[i] => values(g)[i])
        unsafe_push!(f, k, v)
    end
    MA.operate!(canonical, f)
    return f
end

function MA.operate!(
    ::UnsafeAddMul{typeof(*)},
    res::AlgebraElement,
    A::AlgebraElement,
    B::AlgebraElement,
    α = true,
)
    mstr = mstructure(res)
    op = UnsafeAddMul(mstr)
    MA.operate!(op, coeffs(res), coeffs(A), coeffs(B), α)
    return res
end

function MA.operate!(add::UnsafeAdd, res::AlgebraElement, A::AlgebraElement)
    MA.operate!(add, coeffs(res), coeffs(A))
    return res
end

# TODO just push to internal vectors once canonical `does` not just
# call `dropzeros!` but also reorders
function unsafe_push!(a::SparseArrays.SparseVector, k, v)
    a[k] = MA.add!!(a[k], v)
    return a
end

function unsafe_push!(a::Vector, k, v)
    a[k] = MA.add!!(a[k], v)
    return a
end
