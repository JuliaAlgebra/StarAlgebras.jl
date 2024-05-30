function _preallocate_output(X::AlgebraElement, a::Number, op)
    T = MA.promote_operation(op, eltype(X), typeof(a))
    return similar(X, T)
end

function _preallocate_output(X::AlgebraElement, Y::AlgebraElement, op)
    T = MA.promote_operation(op, eltype(X), eltype(Y))
    if coeffs(Y) isa DenseArray # what a hack :)
        return similar(Y, T)
    end
    return similar(X, T)
end

# module structure:

Base.:*(X::AlgebraElement, a::Number) = a * X
Base.:(/)(X::AlgebraElement, a::Number) = inv(a) * X
Base.:(//)(X::AlgebraElement, a::Number) = 1 // a * X

function Base.:-(X::AlgebraElement)
    return MA.operate_to!(_preallocate_output(X, -1, *), -, X)
end
function Base.:*(a::Number, X::AlgebraElement)
    return MA.operate_to!(_preallocate_output(X, a, *), *, X, a)
end
function Base.:div(X::AlgebraElement, a::Number)
    return MA.operate_to!(_preallocate_output(X, a, div), div, X, a)
end

function Base.:+(X::AlgebraElement, Y::AlgebraElement)
    return MA.operate_to!(_preallocate_output(X, Y, +), +, X, Y)
end
function Base.:-(X::AlgebraElement, Y::AlgebraElement)
    return MA.operate_to!(_preallocate_output(X, Y, -), -, X, Y)
end
function Base.:*(X::AlgebraElement, Y::AlgebraElement)
    return MA.operate_to!(_preallocate_output(X, Y, *), *, X, Y)
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
    @assert parent(res) === parent(X)
    @assert parent(X) === parent(Y)
    MA.operate_to!(coeffs(res), +, coeffs(X), coeffs(Y))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(-),
    X::AlgebraElement,
    Y::AlgebraElement,
)
    @assert parent(res) == parent(X) == parent(Y)
    MA.operate_to!(coeffs(res), -, coeffs(X), coeffs(Y))
    return res
end

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(*),
    X::AlgebraElement,
    Y::AlgebraElement,
)
    @assert parent(res) == parent(X) == parent(Y)
    mstr = mstructure(basis(res))
    MA.operate_to!(coeffs(res), mstr, coeffs(X), coeffs(Y))
    return res
end

# TODO just push to internal vectors once canonical is implemented for SparseVector
function unsafe_push!(a::SparseArrays.SparseVector, k, v)
    a[k] = MA.add!!(a[k], v)
    return v
end

function unsafe_push!(a::AlgebraElement, k, v)
    unsafe_push!(coeffs(a), basis(a)[k], v)
    return a
end
