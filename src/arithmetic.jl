# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

_coeff_type(::Type{A}) where {A<:AlgebraElement} = eltype(A)
_coeff_type(a::Type) = a

function _product_coefficient_type(
    ::Type{<:AlgebraElement{T,A}},
    ::Type{Y},
) where {T,A,Y}
    M = MA.promote_operation(mstructure, A)
    C = MA.promote_operation(*, T, _coeff_type(Y))
    return MA.promote_operation(*, value_type(M), C)
end

function algebra_promote_operation(op, args::Vararg{Type,N}) where {N}
    if args[1] <: AlgebraElement &&
       args[2] <: AlgebraElement &&
       MA.promote_operation(parent, args[1]) !==
       MA.promote_operation(parent, args[2])
        return MA.promote_operation(op, _promoted_types(args...)...)
    end
    T = if op === (*) && args[1] <: AlgebraElement && args[2] <: AlgebraElement
        C = _product_coefficient_type(args...)
        MA.promote_operation(+, C, C)
    else
        MA.promote_operation(op, _coeff_type.(args)...)
    end
    if args[2] <: AlgebraElement &&
       MA.promote_operation(coeffs, args[2]) <: DenseArray # what a hack :)
        return similar_type(args[2], T)
    end
    return similar_type(args[1], T)
end

function MA.promote_operation(
    ::typeof(parent),
    ::Type{<:AlgebraElement{T,A}},
) where {T,A}
    return A
end

# module structure:

function Base.:-(X::AlgebraElement)
    return MA.operate_to!(similar(X, MA.promote_operation(-, eltype(X))), -, X)
end
function _scalar_lmul(a, X::AlgebraElement)
    P = MA.promote_operation(*, typeof(a), typeof(X))
    return MA.operate_to!(MA.zero_for(P, (X,)), *, a, X)
end
function _scalar_right(op, X::AlgebraElement, a)
    R = eltype(MA.promote_operation(op, typeof(X), typeof(a)))
    return MA.operate_to!(similar(X, R), op, X, a)
end
# Algebra-valued coefficients are scalars for their containing algebra.
# These signatures resolve the overlap with algebra-by-algebra products.
function Base.:*(a::T, X::AlgebraElement{T}) where {T<:AlgebraElement}
    return _scalar_lmul(a, X)
end
function Base.:*(X::AlgebraElement{T}, a::T) where {T<:AlgebraElement}
    return _scalar_right(*, X, a)
end
function MA.promote_operation(
    ::typeof(*),
    ::Type{T},
    ::Type{A},
) where {T<:AlgebraElement,A<:AlgebraElement{T}}
    return similar_type(A, MA.promote_operation(*, T, T))
end
function MA.promote_operation(
    ::typeof(*),
    ::Type{A},
    ::Type{T},
) where {T<:AlgebraElement,A<:AlgebraElement{T}}
    return similar_type(A, MA.promote_operation(*, T, T))
end
# Array-valued coefficients are scalars for their containing algebra, too.
# Resolve the overlap with MA's array-scalar multiplication promotion.
function MA.promote_operation(
    ::typeof(*),
    ::Type{T},
    ::Type{A},
) where {S,T<:AbstractArray{S},A<:AlgebraElement{T}}
    return similar_type(A, MA.promote_operation(*, T, T))
end
function MA.promote_operation(
    ::typeof(*),
    ::Type{A},
    ::Type{T},
) where {S,T<:AbstractArray{S},A<:AlgebraElement{T}}
    return similar_type(A, MA.promote_operation(*, T, T))
end
function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(*),
    a::T,
    X::AlgebraElement{T},
) where {T<:AlgebraElement}
    return map_coefficients_to!(res, Base.Fix1(*, a), X)
end
function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(*),
    X::AlgebraElement{T},
    a::T,
) where {T<:AlgebraElement}
    return map_coefficients_to!(res, Base.Fix2(*, a), X)
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

_promote_operand!(op, output, x) = x

# Replace the destination's parent and keys together, without changing a parent
# shared with other elements. Type-changing promotion belongs to operate!!.
function _promote_operand!(op, output, x::Union{AlgebraElement,Term})
    parent(x) == parent(output) && return x
    promoted, x = promote_bases(output, x)
    if typeof(promoted) !== typeof(output)
        throw(
            ArgumentError(
                "cannot `$op` in place with a different destination type; use `MutableArithmetics.operate!!`",
            ),
        )
    end
    output.coeffs = coeffs(promoted)
    output.parent = parent(promoted)
    return x
end

_promote_operands!(op, output) = ()

function _promote_operands!(op, output, x, xs::Vararg{Any,N}) where {N}
    x = _promote_operand!(op, output, x)
    ys = _promote_operands!(op, output, xs...)
    # Later operands may have enlarged the common basis again.
    return (_promote_operand!(op, output, x), ys...)
end

"""
    map_coefficients(f, X::AlgebraElement; nonzero=false)

Return an algebra element with the same parent as `X` and coefficients mapped
by `f`. The result's coefficient type is inferred from `f` and `eltype(X)`.
Implicit zeros are left zero. Zero results are removed from sparse storage;
set `nonzero=true` when the mapping is known to preserve nonzero coefficients
to skip this removal.

See also [`map_coefficients!`](@ref) and [`map_coefficients_to!`](@ref).
"""
function map_coefficients(f::F, X::AlgebraElement; nonzero = false) where {F}
    output = similar(X, Base.promote_op(f, eltype(X)))
    return map_coefficients_to!(output, f, X; nonzero)
end

# mutable API

function map_coefficients_to!(
    res::AlgebraElement,
    f::F,
    X::AlgebraElement;
    nonzero = false,
) where {F}
    X = _promote_operand!(map_coefficients_to!, res, X)
    map_coefficients_to!(coeffs(res), f, coeffs(X); nonzero)
    return res
end

function map_coefficients_to!(
    res::AlgebraElement,
    f::F,
    t::Term;
    nonzero = false,
) where {F}
    t = _promote_operand!(map_coefficients_to!, res, t)
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

function MA.operate!(::typeof(one), a::AlgebraElement)
    index = basis(a)[one(object(parent(a)))]
    value = one(eltype(a))
    MA.operate!(zero, a)
    coeffs(a)[index] = value
    return a
end

function MA.operate_to!(
    res::AlgebraElement,
    mul::Union{typeof(__lmul),typeof(__rmul)},
    a,
    X::AlgebraElement,
)
    X = _promote_operand!(mul, res, X)
    MA.operate_to!(coeffs(res), mul, a, coeffs(X))
    return res
end

function MA.operate_to!(res::AlgebraElement, ::typeof(-), X::AlgebraElement)
    X = _promote_operand!(-, res, X)
    MA.operate_to!(coeffs(res), -, coeffs(X))
    return res
end

function MA.operate!(::typeof(+), X::AlgebraElement, Y::AlgebraElement)
    return MA.operate_to!(X, +, X, Y)
end

function MA.operate!(::typeof(-), X::AlgebraElement, Y::AlgebraElement)
    return _scalar_add_mul!(MA.sub_mul, X, Y, identity)
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(+),
    X::AlgebraElement,
    Y::AlgebraElement,
)
    X, Y = _promote_operands!(+, res, X, Y)
    MA.operate_to!(coeffs(res), +, coeffs(X), coeffs(Y))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(-),
    X::AlgebraElement,
    Y::AlgebraElement,
)
    X, Y = _promote_operands!(-, res, X, Y)
    MA.operate_to!(coeffs(res), -, coeffs(X), coeffs(Y))
    return res
end

for (L, R) in (
    (AlgebraElement, AlgebraElement),
    (Term, AlgebraElement),
    (AlgebraElement, Term),
)
    @eval function MA.operate_to!(
        res::AlgebraElement,
        ::typeof(*),
        A::$L,
        B::$R,
    )
        A, B = _promote_operands!(*, res, A, B)
        MA.operate_to!(
            coeffs(res),
            mstructure(res),
            _term_coefficients(A),
            _term_coefficients(B),
        )
        return res
    end
end

function MA.operate!(::typeof(*), A::AlgebraElement, B::AlgebraElement)
    B = _promote_operand!(*, A, B)
    # Preserve the factors while the product kernel clears the destination.
    original = copy(A)
    B = coeffs(A) === coeffs(B) ? original : B
    return MA.operate_to!(A, *, original, B)
end

function MA.operate!(
    op::MA.AddSubMul,
    f::AlgebraElement,
    a::AlgebraElement,
    b::AlgebraElement,
)
    a, b = _promote_operands!(op, f, a, b)
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
    f, a, b = _promote_operands!(op, output, f, a, b)
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
function MA.operate!(
    op::MA.AddSubMul,
    f::AlgebraElement,
    a::MA.Scaling,
    b::MA.Scaling,
    args::Vararg{MA.Scaling,N},
) where {N}
    product = *(map(MA.scaling_to_number, (a, b, args...))...)
    return MA.operate!(MA.add_sub_op(op), f, product)
end

function MA.operate!(
    op::MA.AddSubMul,
    f::AlgebraElement,
    a,
    b,
    args::Vararg{Any,N},
) where {N}
    return MA.operate!(MA.add_sub_op(op), f, MA.operate(*, a, b, args...))
end

function MA.operate_to!(
    output::AlgebraElement,
    op::MA.AddSubMul,
    f::AlgebraElement,
    a,
    b,
    args::Vararg{Any,N},
) where {N}
    return MA.operate_to!(
        output,
        MA.add_sub_op(op),
        f,
        MA.operate(*, a, b, args...),
    )
end

function _scalar_add_mul!(
    op,
    f::AlgebraElement,
    g::AlgebraElement,
    scale::F,
) where {F}
    g = _promote_operand!(op, f, g)
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
