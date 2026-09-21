# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

function MA.promote_operation(
    ::typeof(+),
    ::Type{F},
    ::Type{T},
) where {F<:AlgebraElement,T<:Term}
    return similar_type(F, MA.promote_operation(+, eltype(F), _coeff_type(T)))
end

function MA.operate!(::typeof(+), f::AlgebraElement, t::Term)
    _assert_same_basis(+, f, t)
    iszero(t) && return f
    _add_term!(coeffs(f), t)
    return f
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

function MA.promote_operation(
    op::MA.AddSubMul,
    ::Type{F},
    ::Type{A},
    ::Type{B},
) where {F<:AlgebraElement,A<:Term,B<:Term}
    T = MA.promote_operation(op, eltype(F), _coeff_type(A), _coeff_type(B))
    return similar_type(F, T)
end

function MA.operate!(op::MA.AddSubMul, f::AlgebraElement, a::Term, b::Term)
    c = coeffs(f)
    lt = c isa SparseCoefficients ? c.isless : isless
    # A read-only single entry reuses the term-by-element product kernels.
    single = SparseCoefficients((b.index,), (coefficient(b),), lt)
    return MA.operate!(op, f, a, AlgebraElement(single, parent(b)))
end

# Keep coefficient and basis multiplication in operand order on both sides.
for (L, R, left) in
    ((:Term, :AlgebraElement, true), (:AlgebraElement, :Term, false))
    G = left ? :R : :L
    t, g = left ? (:a, :b) : (:b, :a)
    @eval begin
        function MA.promote_operation(
            op::MA.AddSubMul,
            ::Type{F},
            ::Type{L},
            ::Type{R},
        ) where {F<:AlgebraElement,L<:$L,R<:$R}
            T = MA.promote_operation(*, _coeff_type(L), _coeff_type(R))
            return MA.promote_operation(
                MA.add_sub_op(op),
                F,
                similar_type($G, T),
            )
        end

        function MA.operate!(op::MA.AddSubMul, f::AlgebraElement, a::$L, b::$R)
            return _term_add_mul!(op, f, $t, $g, Val($left))
        end
    end
end

function MA.operate_to!(
    output::AlgebraElement,
    op::MA.AddSubMul,
    f::AlgebraElement,
    a::Union{AlgebraElement,Term},
    b::Union{AlgebraElement,Term},
)
    _assert_same_basis(op, output, f)
    _assert_same_basis(op, output, a)
    _assert_same_basis(op, output, b)
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
        MA.operate!(UnsafeAdd(), output, f)
        MA.operate!(canonical, coeffs(output))
    end
    return MA.operate!(op, output, a, b)
end

_term_product_style(ms, f, g) = GeneralTermProduct()
function _term_product_style(ms, f::SparseCoefficients, g::SparseCoefficients)
    if f.isless == g.isless
        return term_product_style(ms, f.isless)
    end
    return GeneralTermProduct()
end

function _term_add_mul!(op::MA.AddSubMul, f::AlgebraElement, t::Term, g, left)
    _assert_same_basis(op, f, t)
    _assert_same_basis(op, f, g)
    iszero(t) && return f
    ms = mstructure(f)
    c, d = coeffs(f), coeffs(g)
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
    return _merge_coefficients!(f, f, g, product)
end
