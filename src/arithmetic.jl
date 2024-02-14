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
    a.coeffs .= zero(eltype(coeffs(a)))
    return a
end

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

function mul!(
    res::AbstractVector,
    X::AbstractVector,
    Y::AbstractVector,
    ms::MultiplicativeStructure,
)
    res = (res === X || res === Y) ? zero(res) : (res .= zero(eltype(res)))
    return fmac!(res, X, Y, ms)
end

function mul!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    res = (res === X || res === Y) ? zero(res) : zero!(res)
    return fmac!(res, X, Y)
end

function fmac!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    @assert parent(res) === parent(X) === parent(Y)
    fmac!(coeffs(res), coeffs(X), coeffs(Y), parent(res).mstructure)
    return res
end

_nzpairs(v::AbstractVector) = pairs(v)
_nzpairs(v::AbstractSparseVector) =
    zip(SparseArrays.nonzeroinds(v), SparseArrays.nonzeros(v))

function fmac!(
    res::AbstractVector,
    X::AbstractVector,
    Y::AbstractVector,
    mstr::MultiplicativeStructure,
)
    @assert res !== X
    @assert res !== Y
    lx, ly = size(mstr)
    @assert all(iszero, @view(X[lx+1:end]))
    @assert all(iszero, @view(Y[ly+1:end]))
    for iy in 1:ly
        y = Y[iy]
        iszero(y) && continue
        for ix in 1:lx
            x = X[ix]
            iszero(x) && continue
            res[mstr[ix, iy]] += X[ix] * y
        end
    end
    return res
end

function fmac!(
    res::AbstractVector,
    X::AbstractSparseVector,
    Y::AbstractSparseVector,
    mstr::MultiplicativeStructure,
)
    @assert res !== X
    @assert res !== Y
    for iy in SparseArrays.nonzeroinds(Y)
        y = Y[iy]
        for ix in SparseArrays.nonzeroinds(X)
            res[mstr[ix, iy]] += X[ix] * y
        end
    end
    return res
end
