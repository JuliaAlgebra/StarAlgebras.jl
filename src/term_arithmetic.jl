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
            return _term_add_mul!(op, f, $t, $g, Val($left))
        end

        function MA.operate_to!(
            output::AlgebraElement,
            op::MA.AddSubMul,
            f::AlgebraElement,
            a::$L,
            b::$R,
        )
            a, b = _prepare_fused_output!(output, op, f, a, b)
            return _term_add_mul!(op, output, $t, $g, Val($left))
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

function _term_add_mul!(op::MA.AddSubMul, f::AlgebraElement, t::Term, g, left)
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
    _term_add_mul!(style, op, c, ms, t, d, left)
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
) where {left}
    key_product = left ? Base.Fix1(ms, t.index) : Base.Fix2(ms, t.index)
    value_product =
        left ? Base.Fix1(*, coefficient(t)) : Base.Fix2(*, coefficient(t))
    add = MA.add_sub_op(op)
    product = function (p)
        k, scale = only(nonzero_pairs(key_product(first(p))))
        return k => add(scale * value_product(last(p)))
    end
    return _merge_coefficients!(f, f, g, product, identity)
end
