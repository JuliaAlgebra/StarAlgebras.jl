# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

# Polynomial frontends can extend a finite basis before adding a constant.
_prepare_scalar_sum(p) = p
MA.promote_operation(::typeof(_prepare_scalar_sum), ::Type{P}) where {P} = P

_coeff_type(::Type{<:Term{T}}) where {T} = T

function _scalar_constant(op, p, a, ::Val{left}) where {left}
    T = _coeff_type(typeof(p))
    C =
        left ? MA.promote_operation(op, typeof(a), T) :
        MA.promote_operation(op, T, typeof(a))
    # A scalar such as a variable index needs a coefficient type before zero
    # initialization or negation. Use the same promotion for every scalar family.
    return _constant_term(p, convert(C, a))
end

function _scalar_sum(op, p, a, left::Val{L}) where {L}
    p = _prepare_scalar_sum(p)
    t = _scalar_constant(op, p, a, left)
    return L ? op(t, p) : op(p, t)
end

function _scalar_sum_to!(output, op, p, a, left::Val{L}) where {L}
    p = _promote_operand!(op, output, p)
    t = _scalar_constant(op, output, a, left)
    x, y = L ? (t, p) : (p, t)
    MA.operate_to!(
        coeffs(output),
        op,
        _term_coefficients(x),
        _term_coefficients(y),
    )
    return output
end

# Generate the same methods for the coefficient type and registered scalar
# families. The coefficient-type methods resolve intersections between families.
function _scalar_arithmetic(S)
    code = quote
        function Base.:*(a::$S, t::Term{T}) where {T}
            return Term(parent(t), MA.copy_if_mutable(t.index), a * coefficient(t))
        end

        function Base.:*(t::Term{T}, a::$S) where {T}
            return Term(parent(t), MA.copy_if_mutable(t.index), coefficient(t) * a)
        end

        function MA.promote_operation(
            ::typeof(*),
            ::Type{C},
            ::Type{Term{T,A,I}},
        ) where {T,A,I,C<:$S}
            return Term{MA.promote_operation(*, C, T),A,I}
        end

        function MA.promote_operation(
            ::typeof(*),
            ::Type{Term{T,A,I}},
            ::Type{C},
        ) where {T,A,I,C<:$S}
            return Term{MA.promote_operation(*, T, C),A,I}
        end

        Base.:*(a::$S, p::AlgebraElement{T}) where {T} = _scalar_lmul(a, p)

        function MA.promote_operation(
            ::typeof(*),
            ::Type{C},
            ::Type{P},
        ) where {T,P<:AlgebraElement{T},C<:$S}
            return similar_type(P, MA.promote_operation(*, C, T))
        end

        function MA.promote_operation(
            op::Union{typeof(*),typeof(/),typeof(//),typeof(div)},
            ::Type{P},
            ::Type{C},
        ) where {T,P<:AlgebraElement{T},C<:$S}
            R =
                op === (//) ? Base.promote_op(op, T, C) :
                MA.promote_operation(op, T, C)
            return similar_type(P, R)
        end

        function MA.operate_to!(
            output::AlgebraElement,
            ::typeof(*),
            a::$S,
            p::AlgebraElement{T},
        ) where {T}
            return map_coefficients_to!(output, Base.Fix1(*, a), p)
        end

        function MA.operate_to!(
            output::AlgebraElement,
            op::Union{typeof(*),typeof(/),typeof(//),typeof(div)},
            p::AlgebraElement{T},
            a::$S,
        ) where {T}
            return map_coefficients_to!(output, Base.Fix2(op, a), p)
        end

        function MA.operate!(
            op::Union{typeof(*),typeof(/),typeof(//),typeof(div)},
            p::AlgebraElement{T},
            a::$S,
        ) where {T}
            return MA.operate_to!(p, op, p, a)
        end

        function Base.:*(a::$S, c::AbstractCoefficients{K,T}) where {K,T}
            output = similar(c, MA.promote_operation(*, typeof(a), T))
            return MA.operate_to!(output, *, a, c)
        end

        function MA.operate_to!(
            output::AbstractCoefficients,
            ::typeof(*),
            a::$S,
            c::AbstractCoefficients{K,T},
        ) where {K,T}
            return map_coefficients_to!(output, Base.Fix1(*, a), c)
        end

        function MA.operate_to!(
            output::AbstractCoefficients,
            op::Union{typeof(*),typeof(/),typeof(//),typeof(div)},
            c::AbstractCoefficients{K,T},
            a::$S,
        ) where {K,T}
            return map_coefficients_to!(output, Base.Fix2(op, a), c)
        end

        function MA.operate!(::typeof(+), p::AlgebraElement{T}, a::$S) where {T}
            return MA.operate!(+, p, _constant_term(p, a))
        end

        function MA.operate!(::typeof(-), p::AlgebraElement{T}, a::$S) where {T}
            return MA.operate!(+, p, -convert(T, a))
        end
    end
    for op in (:*, :/, ://, :div)
        push!(
            code.args,
            quote
                function Base.$op(p::AlgebraElement{T}, a::$S) where {T}
                    return _scalar_right($op, p, a)
                end

                function Base.$op(
                    c::AbstractCoefficients{K,T},
                    a::$S,
                ) where {K,T}
                    R =
                        $op === (//) ? Base.promote_op($op, T, typeof(a)) :
                        MA.promote_operation($op, T, typeof(a))
                    return MA.operate_to!(similar(c, R), $op, c, a)
                end
            end,
        )
    end
    for left in (true, false)
        L, R =
            left ? (S, :(Union{AlgebraElement{T},Term{T}})) :
            (:(Union{AlgebraElement{T},Term{T}}), S)
        scalar, p = left ? (:a, :b) : (:b, :a)
        fix = left ? Base.Fix1 : Base.Fix2
        push!(
            code.args,
            quote
                function MA.operate!(
                    op::MA.AddSubMul,
                    output::AlgebraElement,
                    a::$L,
                    b::$R,
                ) where {T}
                    return _scalar_add_mul!(op, output, $p, $fix(*, $scalar))
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
            end,
        )
        for op in (:+, :-)
            push!(
                code.args,
                quote
                    function Base.$op(a::$L, b::$R) where {T}
                        return _scalar_sum($op, $p, $scalar, Val($left))
                    end

                    function MA.promote_operation(
                        ::typeof($op),
                        ::Type{A},
                        ::Type{B},
                    ) where {T,A<:$L,B<:$R}
                        P, C = $left ? (B, A) : (A, B)
                        P = MA.promote_operation(_prepare_scalar_sum, P)
                        R =
                            $left ? MA.promote_operation($op, C, T) :
                            MA.promote_operation($op, T, C)
                        return similar_type(_algebra_element_type(P), R)
                    end

                    function MA.operate_to!(
                        output::AlgebraElement,
                        ::typeof($op),
                        a::$L,
                        b::$R,
                    ) where {T}
                        return _scalar_sum_to!(
                            output,
                            $op,
                            $p,
                            $scalar,
                            Val($left),
                        )
                    end
                end,
            )
        end
    end
    return code
end

"""
    @scalar_arithmetic Scalar

Allow values of `Scalar` to act as coefficients of [`Term`](@ref) and
[`AlgebraElement`](@ref), independently of their current coefficient type.
StarAlgebras uses this registration for its built-in scalar arithmetic too.

This defines addition, subtraction and multiplication on both sides; right
division (`/`, `//`, `div`) of algebra elements; scalar construction via an
algebra; and the corresponding MutableArithmetics promotion and accumulation
methods. Coefficient containers use the same scalar multiplication and division.
Each operation requires the corresponding coefficient operation and its
`MutableArithmetics.promote_operation` to be defined. Multiplication preserves
operand order. Addition and subtraction require a basis containing the identity.

Use this macro at module scope, typically in a package extension. Register
disjoint scalar families; values of the existing coefficient type are already
supported without registration. This does not make arbitrary objects scalars.
"""
macro scalar_arithmetic(Scalar)
    return quote
        $(_scalar_arithmetic(:(Union{T,$(esc(Scalar))})))
        function (A::AbstractStarAlgebra)(x::$(esc(Scalar)))
            return __coerce(A, (one(object(A)) => x))
        end
    end
end

@eval $(_scalar_arithmetic(:T))
@scalar_arithmetic Number
