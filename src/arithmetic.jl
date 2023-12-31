# module structure:

function Base.:*(a::Number, X::AlgebraElement)
    T = Base._return_type(*, Tuple{eltype(X),typeof(a)})
    return mul!(similar(X, T), X, a)
end

Base.:*(X::AlgebraElement, a::Number) = a * X
Base.:(/)(X::AlgebraElement, a::Number) = inv(a) * X

# TODO: handle this through mul!?
Base.:(//)(X::AlgebraElement, a::Number) = AlgebraElement(coeffs(X) .// a, parent(X))
Base.:div(X::AlgebraElement, a::Number) = AlgebraElement(div.(coeffs(X), a), parent(X))

# ring structure:
Base.:-(X::AlgebraElement) = neg!(similar(X), X)

function _preallocate_output(X::AlgebraElement, Y::AlgebraElement, op)
    T = Base._return_type(op, Tuple{eltype(X),eltype(Y)})
    if coeffs(Y) isa DenseArray
        return similar(Y, T)
    end
    return similar(X, T)
end

Base.:+(X::AlgebraElement, Y::AlgebraElement) = add!(_preallocate_output(X, Y, +), X, Y)
Base.:-(X::AlgebraElement, Y::AlgebraElement) = sub!(_preallocate_output(X, Y, -), X, Y)
Base.:*(X::AlgebraElement, Y::AlgebraElement) = mul!(_preallocate_output(X, Y, *), X, Y)

Base.:^(a::AlgebraElement, p::Integer) = Base.power_by_squaring(a, p)

# mutable API; TODO: replace with MutableArithmetic

function zero!(a::AlgebraElement)
    MA.operate!(zero, a.coeffs)
    return a
end

MA.operate!(::typeof(zero), v::SparseVector) = (v .*= 0; v)

function neg!(res::AlgebraElement, X::AlgebraElement)
    @assert parent(res) === parent(X)
    res.coeffs .= .-coeffs(X)
    return res
end

function add!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    @assert parent(res) === parent(X)
    @assert parent(X) === parent(Y)
    if res === X
        for (idx, y) in _nzpairs(coeffs(Y))
            res[idx] += y
        end
    elseif res === Y
        for (idx, x) in _nzpairs(coeffs(X))
            res[idx] += x
        end
    else
        zero!(res)
        for (idx, x) in _nzpairs(coeffs(X))
            res[idx] += x
        end
        for (idx, y) in _nzpairs(coeffs(Y))
            res[idx] += y
        end
    end
    return res
end

function sub!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    @assert parent(res) === parent(X) === parent(Y)
    neg!(res, Y)
    add!(res, res, X)
    return res
end

function mul!(res::AlgebraElement, X::AlgebraElement, a::Number)
    @assert parent(res) === parent(X)
    if res !== X
        zero!(res)
    end
    for (idx, x) in _nzpairs(coeffs(X))
        res.coeffs[idx] = a * x
    end
    return res
end

function mul!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    @assert parent(res) === parent(X) === parent(Y)
    mstr = mstructure(basis(parent(res)))
    mul!(mstr, coeffs(res), coeffs(X), coeffs(Y))
    return res
end

_nzpairs(v::AbstractVector) = pairs(v)
_nzpairs(v::AbstractSparseVector) =
    zip(SparseArrays.nonzeroinds(v), SparseArrays.nonzeros(v))
