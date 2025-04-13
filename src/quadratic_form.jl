# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

"""
    QuadraticForm(Q)
    QuadraticForm{star}(Q)
A simple wrapper for representing a quadratic form.

`QuadraticForm(Q)` represents an honest quadratic form, however in the context
of `*`-algebras a more canonical object is a **sesquilinear form** which can be
constructed as `QuadraticForm{star}(Q)`.

Both objects are defined by matrix coefficients with respect to their bases `b`,
i.e. their values at `i`-th and `j`-th basis elements:
    `star.(b)·Q·b = ∑ᵢⱼ Q[i,j]·star(b[i])·b[j]`.

# Necessary methods:
 * `basis(Q)` - a **finite** basis `b` of the quadratic form;
 * `Base.getindex(Q, i::T, j::T)` - the value at `i`-th and `j`-th basis elements
   where `T` is the type of indicies of `basis(Q)`;
 * `Base.eltype(Q)` - the type of `Q[i,j]·star(b[i])·b[j]`.
"""
struct QuadraticForm{T,involution}
    Q::T
    QuadraticForm(Q) = new{typeof(Q),identity}(Q)
    QuadraticForm{star}(Q) = new{typeof(Q),star}(Q)
end

Base.eltype(qf::QuadraticForm) = eltype(qf.Q)
basis(qf::QuadraticForm) = basis(qf.Q)
Base.getindex(qf::QuadraticForm, i::T, j::T) where {T} = qf.Q[i, j]

function MA.operate_to!(
    res::AlgebraElement,
    ::typeof(+),
    Q::QuadraticForm,
)
    MA.operate!(zero, res)
    MA.operate!(UnsafeAdd(), res, Q)
    MA.operate!(canonical, res)
    return res
end

function MA.operate!(
    ::UnsafeAdd,
    res::AlgebraElement,
    Q::QuadraticForm{T,ε},
) where {T,ε}
    op = UnsafeAddMul(mstructure(basis(res)))
    for (i, b1) in pairs(basis(Q))
        b1★ = ε(b1)
        for (j, b2) in pairs(basis(Q))
            MA.operate!(op, res, b1★, b2, Q[i, j])
        end
    end
end

function MA.operate!(op::UnsafeAddMul, res::AlgebraElement{A,T}, b1::T, b2::T, a) where {A, T}
    b = basis(res)
    return MA.operate!(op, coeffs(res), mstructure(b)(b[b1], b[b2]), a)
end

"""
    AlgebraElement(qf::QuadraticForm, A::AbstractStarAlgebra)
Construct an algebra element in `A` representing quadratic form `qf`.

!!! warning
    It is assumed that all basis elements of `qf` belong to `A`, or at least
    that `keys` of their coefficients can be found in the basis of `A`.
"""
function AlgebraElement(qf::QuadraticForm, A::AbstractStarAlgebra)
    @assert all(b -> parent(b) == A, basis(qf))
    res = zero(eltype(qf), A)
    MA.operate_to!(res, +, qf)
    return res
end

(A::AbstractStarAlgebra)(qf::QuadraticForm) = AlgebraElement(qf, A)
