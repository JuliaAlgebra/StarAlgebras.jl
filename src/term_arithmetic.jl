# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

function MA.operate!(::typeof(+), f::AlgebraElement, t::Term)
    t = _promote_operand!(+, f, t)
    iszero(t) && return f
    _add_term!(coeffs(f), t)
    return f
end

function MA.operate!(::typeof(-), f::AlgebraElement, t::Term)
    return MA.operate!(+, f, -t)
end

function MA.operate_to!(
    output::AlgebraElement,
    op::Union{typeof(+),typeof(-)},
    f::AlgebraElement,
    t::Term,
)
    fused = op === (+) ? MA.add_mul : MA.sub_mul
    return MA.operate_to!(output, fused, f, true, t)
end

function _add_term!(c, t::Term)
    c[t.index] += coefficient(t)
    MA.operate!(canonical, c)
    return c
end

function _add_term!(c::SparseCoefficients, t::Term)
    # Reuse the sorted merge without allocating a single-term polynomial.
    # Inserted mutable keys and coefficients must remain independent of `t`.
    single = SparseCoefficients(
        (MA.copy_if_mutable(t.index),),
        (MA.copy_if_mutable(coefficient(t)),),
        c.isless,
    )
    return MA.operate!(+, c, single)
end

# Keep coefficient and basis multiplication in operand order on both sides.
for (L, R, left) in (
    (:Term, :AlgebraElement, true),
    (:AlgebraElement, :Term, false),
    (:Term, :Term, true),
)
    t, g = left ? (:a, :b) : (:b, :a)
    @eval begin
        function MA.operate!(op::MA.AddSubMul, f::AlgebraElement, a::$L, b::$R)
            buffer = MA.buffer_for(op, typeof(f), typeof(a), typeof(b))
            return MA.buffered_operate!(buffer, op, f, a, b)
        end

        function MA.buffer_for(
            op::MA.AddSubMul,
            ::Type{F},
            ::Type{A},
            ::Type{B},
        ) where {F<:AlgebraElement,A<:$L,B<:$R}
            return MA.buffer_for(op, _coefficient_type.((F, A, B))...)
        end

        function MA.buffered_operate!(
            buffer,
            op::MA.AddSubMul,
            f::AlgebraElement,
            a::$L,
            b::$R,
        )
            return _term_add_mul!(op, f, $t, $g, Val($left), buffer)
        end

        function MA.operate_to!(
            output::AlgebraElement,
            op::MA.AddSubMul,
            f::AlgebraElement,
            a::$L,
            b::$R,
        )
            buffer = MA.buffer_for(op, typeof(f), typeof(a), typeof(b))
            return MA.buffered_operate_to!(buffer, output, op, f, a, b)
        end

        function MA.buffered_operate_to!(
            buffer,
            output::AlgebraElement,
            op::MA.AddSubMul,
            f::AlgebraElement,
            a::$L,
            b::$R,
        )
            a, b = _prepare_fused_output!(output, op, f, a, b)
            return _term_add_mul!(op, output, $t, $g, Val($left), buffer)
        end
    end
end

_term_product_style(ms, f, g) = GeneralTermProduct()
function _term_product_style(ms, f::SparseCoefficients, g::SparseCoefficients)
    if f.isless == g.isless
        return term_product_style(ms, f.isless)
    end
    return GeneralTermProduct()
end

_coefficient_type(::Type{<:Union{AlgebraElement{T},Term{T}}}) where {T} = T

# `c` multiplies every term so it must not be one of the coefficients of `f`
# that are modified in place.
function _not_in_values(c, f)
    if MA.mutability(typeof(c)) isa MA.IsMutable &&
       any(Base.Fix1(===, c), values(f))
        return MA.copy_if_mutable(c)
    end
    return c
end

function _term_add_mul!(
    op::MA.AddSubMul,
    f::AlgebraElement,
    t::Term,
    g,
    left,
    buffer = nothing,
)
    t, g = _promote_operands!(op, f, t, g)
    iszero(t) && return f
    ms = mstructure(f)
    c = coeffs(f)
    d =
        g isa Term ?
        SparseCoefficients(
            (g.index,),
            (coefficient(g),),
            comparable(basis(f)),
        ) : coeffs(g)
    style = _term_product_style(ms, c, d)
    _term_add_mul!(style, op, c, ms, t, d, left, buffer)
    return f
end

function _term_add_mul!(
    ::GeneralTermProduct,
    op::MA.AddSubMul,
    f,
    ms,
    t,
    g,
    ::Val{left},
    _,
) where {left}
    # Accumulation may grow or overwrite the storage being read.
    g = f === g ? MA.mutable_copy(g) : g
    c = SparseCoefficients((t.index,), (MA.copy_if_mutable(coefficient(t)),))
    a, b = left ? (c, g) : (g, c)
    α = op === MA.add_mul ? true : -1
    MA.operate!(UnsafeAddMul(ms), f, a, b, α)
    MA.operate!(canonical, f)
    return f
end

function _term_add_mul!(
    ::OrderedTermProduct,
    op::MA.AddSubMul,
    f::SparseCoefficients,
    ms,
    t,
    g,
    ::Val{left},
    buffer,
) where {left}
    key_product = left ? Base.Fix1(ms, t.index) : Base.Fix2(ms, t.index)
    c = _not_in_values(coefficient(t), f)
    factors = left ? Base.Fix1(tuple, c) : Base.Fix2(tuple, c)
    add = MA.add_sub_op(op)
    # Keep the coefficient of `g` so that the product is only computed once
    # we know whether it is added to a coefficient of `f`.
    keyed = function (p)
        k, scale = only(nonzero_pairs(key_product(first(p))))
        return k => (scale, last(p))
    end
    # `f` owns its coefficients so they can be modified in place.
    combine = function (a, b)
        scale, v = last(b)
        l, r = factors(v)
        if isone(scale)
            value = MA.buffered_operate!!(buffer, op, last(a), l, r)
        else
            value = MA.operate!!(op, last(a), scale * l, r)
        end
        return first(a) => value
    end
    new_term = function (b)
        scale, v = last(b)
        return first(b) => add(scale * *(factors(v)...))
    end
    _merge_sorted!(
        _CoefficientPairs(f),
        _CoefficientPairs(f),
        _CoefficientPairs(g);
        lt = f.isless,
        by = first,
        map_b = keyed,
        combine,
        filter = p -> !iszero(last(p)),
        only_b = new_term,
    )
    return f
end
