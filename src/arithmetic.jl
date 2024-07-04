_coeff_type(::Type{A}) where {A<:AlgebraElement} = eltype(A)
_coeff_type(a::Type) = a
_coeff_type(a) = _coeff_type(typeof(a))

function algebra_promote_operation(op, args::Vararg{Type,N}) where {N}
    T = MA.promote_operation(op, _coeff_type.(args)...)
    if args[2] <: AlgebraElement && MA.promote_operation(coeffs, args[2]) <: DenseArray # what a hack :)
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

Base.:*(X::AlgebraElement, a::Number) = a * X
Base.:(/)(X::AlgebraElement, a::Number) = inv(a) * X
Base.:(//)(X::AlgebraElement, a::Number) = 1 // a * X

function Base.:-(X::AlgebraElement)
    return MA.operate_to!(_preallocate_output(*, X, -1), -, X)
end
function MA.promote_operation(::typeof(*), ::Type{T}, ::Type{A}) where {T<:Number, A<:AlgebraElement}
    return algebra_promote_operation(*, A, T)
end
function Base.:*(a::Any, X::AlgebraElement)
    return MA.operate_to!(_preallocate_output(*, X, a), *, a, X)
end
function Base.:div(X::AlgebraElement, a::Number)
    return MA.operate_to!(_preallocate_output(div, X, a), div, X, a)
end
function Base.:*(
    a::T,
    X::AlgebraElement{A},
) where {T,O,A<:AbstractStarAlgebra{O,T}}
    return MA.operate_to!(similar(X), Val(:lmul), a, X)
end
function Base.:*(
    X::AlgebraElement{A},
    a::T,
) where {T,O,A<:AbstractStarAlgebra{O,T}}
    return MA.operate_to!(similar(X), Val(:rmul), a, X)
end

for op in [:+, :-, :*]
    @eval begin
        function MA.promote_operation(::typeof($op), ::Type{X}, ::Type{Y}) where {X<:AlgebraElement,Y<:AlgebraElement}
            return algebra_promote_operation($op, X, Y)
        end
        function Base.$op(X::AlgebraElement, Y::AlgebraElement)
            return MA.operate_to!(_preallocate_output($op, X, Y), $op, X, Y)
        end
    end
end

function Base.:*(args::Vararg{AlgebraElement,N}) where {N}
    return MA.operate_to!(_preallocate_output(*, args...), *, args...)
end

Base.:^(a::AlgebraElement, p::Integer) = Base.power_by_squaring(a, p)

# mutable API

function MA.operate!(::typeof(zero), a::AlgebraElement)
    MA.operate!(zero, coeffs(a))
    return a
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(*),
    a::Any,
    X::AlgebraElement,
)
    @assert parent(res) === parent(X)
    MA.operate_to!(coeffs(res), *, a, coeffs(X))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(div),
    X::AlgebraElement,
    a::Number,
)
    @assert parent(res) === parent(X)
    MA.operate_to!(coeffs(res), div, coeffs(X), a)
    return res
end

function MA.operate_to!(res::AlgebraElement, mul::Val, a, X::AlgebraElement)
    @assert parent(res) == parent(X)
    MA.operate_to!(coeffs(res), mul, a, coeffs(X))
    return res
end

function MA.operate_to!(res::AlgebraElement, ::typeof(-), X::AlgebraElement)
    @assert parent(res) === parent(X)
    MA.operate_to!(coeffs(res), -, coeffs(X))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(+),
    X::AlgebraElement,
    Y::AlgebraElement,
)
    @assert parent(res) == parent(X)
    @assert parent(X) == parent(Y)
    MA.operate_to!(coeffs(res), +, coeffs(X), coeffs(Y))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(-),
    X::AlgebraElement,
    Y::AlgebraElement,
)
    @assert parent(res) === parent(X) === parent(Y)
    MA.operate_to!(coeffs(res), -, coeffs(X), coeffs(Y))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(*),
    args::Vararg{AlgebraElement,N},
) where {N}
    for arg in args
        @assert parent(res) == parent(arg)
    end
    mstr = mstructure(basis(res))
    MA.operate_to!(coeffs(res), mstr, coeffs.(args)...)
    return res
end

function MA.operate!(
    ::UnsafeAddMul{typeof(*)},
    res::AlgebraElement,
    args::Vararg{AlgebraElement,N},
) where {N}
    for arg in args
        @assert parent(res) == parent(arg)
    end
    mstr = mstructure(basis(res))
    MA.operate!(UnsafeAddMul(mstr), coeffs(res), coeffs.(args)...)
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
