# module structure:
Base.:*(a::Number, X::AlgebraElement) =
    mul!(similar(X, promote_type(eltype(X), typeof(a))), X, a)
Base.:*(X::AlgebraElement, a::Number) = a * X
Base.:(/)(X::AlgebraElement, a::Number) = inv(a) * X

# TODO: handle this through mul!?
Base.:(//)(X::AlgebraElement, a::Number) = AlgebraElement(coeffs(X) .// a, parent(X))

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

zero!(v::AbstractVector{T}) where {T} = (v .= zero(T); v)

function zero!(a::AlgebraElement)
    a.coeffs .= zero(first(coeffs(a)))
    return a
end

function neg!(res::AlgebraElement, X::AlgebraElement)
    @assert parent(res) === parent(X)
    res.coeffs .= -coeffs(X)
    return res
end

function add!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    @assert parent(res) === parent(X) === parent(Y)
    # res = (res === X || res === Y) ? similar(res) : res
    res.coeffs .= coeffs(X) .+ coeffs(Y)
    return res
end

function sub!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    @assert parent(res) === parent(X) === parent(Y)
    # res = (res === X || res === Y) ? similar(res) : res
    res.coeffs .= coeffs(X) .- coeffs(Y)
    return res
end

function mul!(res::AlgebraElement, X::AlgebraElement, a::Number)
    @assert parent(res) === parent(X)
    # res = (res === X || res === Y) ? similar(res) : res
    res.coeffs .= a .* coeffs(X)
    return res
end

function mul!(
    res::AbstractVector,
    X::AbstractVector,
    Y::AbstractVector,
    ms::MultiplicativeStructure,
)
    res = (res === X || res === Y) ? zero(res) : zero!(res)
    return fmac!(res, X, Y, ms)
end

function mul!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    @assert parent(res) === parent(X) === parent(Y)
    res = (res === X || res === Y) ? zero(res) : zero!(res)
    return fmac!(res, X, Y)
end

function fmac!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    fmac!(coeffs(res), coeffs(X), coeffs(Y), parent(res).mstructure)
    return res
end

_nzidx(v::AbstractVector) = eachindex(v)
_nzidx(v::AbstractSparseVector) = SparseArrays.nonzeroinds(v)

function fmac!(
    res::AbstractVector,
    X::AbstractVector,
    Y::AbstractVector,
    mstr::MultiplicativeStructure,
)
    @inbounds for j in _nzidx(Y)
        Yj = Y[j]
        for i in _nzidx(X)
            res[mstr[i, j]] += X[i] * Yj
        end
    end
    return res
end
