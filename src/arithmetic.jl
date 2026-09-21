# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

_coeff_type(::Type{A}) where {A<:AlgebraElement} = eltype(A)
_coeff_type(::Type{<:Term{T}}) where {T} = T
_coeff_type(a::Type) = a
_coeff_type(a) = _coeff_type(typeof(a))

function algebra_promote_operation(op, args::Vararg{Type,N}) where {N}
    T = MA.promote_operation(op, _coeff_type.(args)...)
    if args[2] <: AlgebraElement &&
       MA.promote_operation(coeffs, args[2]) <: DenseArray # what a hack :)
        return similar_type(args[2], T)
    end
    return similar_type(args[1], T)
end

function _preallocate_output(op, args::Vararg{Any,N}) where {N}
    T = MA.promote_operation(op, _coeff_type.(args)...)
    if args[2] isa AlgebraElement && coeffs(args[2]) isa DenseArray # what a hack :)
        return similar(args[2], T)
    end
    return similar(args[1], T)
end

# module structure:

function Base.:-(X::AlgebraElement)
    return MA.operate_to!(similar(X, MA.promote_operation(-, eltype(X))), -, X)
end
function MA.promote_operation(
    ::typeof(*),
    ::Type{<:Union{T,Number}},
    ::Type{A},
) where {T,A<:AlgebraElement{T}}
    return algebra_promote_operation(*, A, T)
end
function MA.promote_operation(
    op::Union{typeof(*),typeof(/),typeof(//),typeof(div)},
    ::Type{A},
    ::Type{<:Union{T,Number}},
) where {T,A<:AlgebraElement{T}}
    if op === (//)
        return similar_type(A, Base.promote_op(op, T, T))
    end
    return algebra_promote_operation(op, A, T)
end
function _scalar_lmul(a, X::AlgebraElement{T}) where {T}
    c = convert(T, a)
    return MA.operate_to!(_preallocate_output(*, X, c), *, c, X)
end
function _scalar_right(op, X::AlgebraElement{T}, a) where {T}
    c = convert(T, a)
    R = eltype(MA.promote_operation(op, typeof(X), T))
    return MA.operate_to!(similar(X, R), op, X, c)
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
            return MA.operate_to!(_preallocate_output($op, _X, _Y), $op, _X, _Y)
        end
    end
end

Base.:^(a::AlgebraElement, p::Integer) = Base.power_by_squaring(a, p)

# Informative error for the in-place operations below, which require their
# operands to share a basis (they do not promote, unlike `*`, `+`, `-`).
function _assert_same_basis(
    op,
    A::AlgebraElement,
    B::Union{AlgebraElement,Term},
)
    parent(A) == parent(B) && return
    return throw(
        ArgumentError(
            "cannot `$op` algebra elements over different bases in place: their " *
            "bases differ. Bring them to a common basis first, e.g. " *
            "`_A, _B = StarAlgebras.promote_bases(A, B)`, or use the `*`, `+`, `-` " *
            "operators which promote automatically.",
        ),
    )
end

# mutable API

function map_coefficients_to!(
    res::AlgebraElement,
    f::F,
    X::AlgebraElement;
    nonzero = false,
) where {F}
    _assert_same_basis(map_coefficients_to!, res, X)
    map_coefficients_to!(coeffs(res), f, coeffs(X); nonzero)
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
    @assert parent(res) === parent(X)
    MA.operate_to!(coeffs(res), *, convert(T, a), coeffs(X))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    op::Union{typeof(*),typeof(/),typeof(//),typeof(div)},
    X::AlgebraElement{T},
    a::Union{T,Number},
) where {T}
    @assert parent(res) === parent(X)
    MA.operate_to!(coeffs(res), op, coeffs(X), convert(T, a))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    mul::Union{typeof(__lmul),typeof(__rmul)},
    a,
    X::AlgebraElement,
)
    @assert parent(res) == parent(X)
    MA.operate_to!(coeffs(res), mul, a, coeffs(X))
    return res
end

function MA.operate_to!(res::AlgebraElement, ::typeof(-), X::AlgebraElement)
    @assert parent(res) === parent(X)
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
    _assert_same_basis(+, X, Y)
    @assert parent(res) == parent(X)
    MA.operate_to!(coeffs(res), +, coeffs(X), coeffs(Y))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(-),
    X::AlgebraElement,
    Y::AlgebraElement,
)
    _assert_same_basis(-, X, Y)
    @assert parent(res) == parent(X)
    MA.operate_to!(coeffs(res), -, coeffs(X), coeffs(Y))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(*),
    A::AlgebraElement,
    B::AlgebraElement,
)
    _assert_same_basis(*, A, B)
    @assert parent(res) == parent(A)
    mstr = mstructure(res)
    MA.operate_to!(coeffs(res), mstr, coeffs(A), coeffs(B), true)
    return res
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
