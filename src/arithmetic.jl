_coeff_type(X::AlgebraElement) = eltype(X)
_coeff_type(a) = typeof(a)

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
function Base.:*(a::Number, X::AlgebraElement)
    return MA.operate_to!(_preallocate_output(*, X, a), *, X, a)
end
function Base.:div(X::AlgebraElement, a::Number)
    return MA.operate_to!(_preallocate_output(div, X, a), div, X, a)
end

function Base.:+(X::AlgebraElement, Y::AlgebraElement)
    return MA.operate_to!(_preallocate_output(+, X, Y), +, X, Y)
end
function Base.:-(X::AlgebraElement, Y::AlgebraElement)
    return MA.operate_to!(_preallocate_output(-, X, Y), -, X, Y)
end
function Base.:*(X::AlgebraElement, Y::AlgebraElement)
    return MA.operate_to!(_preallocate_output(*, X, Y), *, X, Y)
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
    X::AlgebraElement,
    a::Number,
)
    @assert parent(res) === parent(X)
    MA.operate_to!(coeffs(res), *, coeffs(X), a)
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
        if arg isa AlgebraElement
            @assert parent(res) == parent(arg)
        end
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
        if arg isa AlgebraElement
            @assert parent(res) == parent(arg)
        end
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
