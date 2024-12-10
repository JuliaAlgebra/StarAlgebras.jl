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
    res,
    ms::MultiplicativeStructure,
    Q::QuadraticForm{T,ε},
) where {T,ε}
    MA.operate!(zero, res)
    op = UnsafeAddMul(ms)
    for (i, b1) in pairs(basis(Q))
        b1★ = ε(b1)
        for (j, b2) in pairs(basis(Q))
            MA.operate!(op, res, coeffs(b1★), coeffs(b2), Q[i, j])
        end
    end
    MA.operate!(canonical, res)
    return res
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
    MA.operate_to!(coeffs(res), mstructure(A), qf)
    return res
end

(A::AbstractStarAlgebra)(qf::QuadraticForm) = AlgebraElement(qf, A)
