function Base.show(io::IO, A::AbstractStarAlgebra)
    ioc = IOContext(io, :limit => true, :compact => true)
    return print(ioc, "*-algebra of ", object(A))
end

__prints_with_minus(::Any) = false
__prints_with_minus(x::Real) = x < 0
__needs_parens(::Any) = false
__needs_parens(a::AlgebraElement) = true

function _coeff_elt_print(io, c, elt)
    print(io, c, '·')
    __needs_parens(elt) && print(io, '(')
    print(io, elt)
    __needs_parens(elt) && print(io, ')')
    return
end

function Base.show(io::IO, a::AlgebraElement)
    A = parent(a)
    if iszero(a)
        T = valtype(coeffs(a))
        _coeff_elt_print(io, zero(T), first(basis(A)))
    else
        first = true
        for (idx, value) in pairs(coeffs(a))
            iszero(value) && continue
            c, elt = value, basis(A)[idx]
            if first
                _coeff_elt_print(io, c, elt)
                first = false
            else
                if __prints_with_minus(c)
                    print(io, ' ')
                else
                    print(io, ' ', '+')
                end
                _coeff_elt_print(io, c, elt)
            end
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", mstr::DiracMStructure)
    print(io, "DiracMStructure of ", object(basis(mstr)), " over ")
    print(io, basis(mstr) isa ImplicitBasis ? "implicit" : "explicit")
    print(io, " basis with ")

    print(io, Base.haslength(basis(mstr)) ? length(basis(mstr)) : "indefinite number of")
    print(io, " elements")
end
