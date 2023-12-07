function Base.show(io::IO, A::AbstractStarAlgebra)
    ioc = IOContext(io, :limit => true, :compact => true)
    return print(ioc, "*-algebra of ", object(A))
end

__prints_with_minus(::Any) = false
__prints_with_minus(x::Real) = x < 0
__needs_parens(::Any) = false
__needs_parens(a::AlgebraElement) = true

function _coeff_elt_print(io, c, elt)
    print(io, c, "·")
    __needs_parens(elt) && print(io, '(')
    print(io, elt)
    __needs_parens(elt) && print(io, ')')
    return
end

function Base.show(io::IO, a::AlgebraElement)
    A = parent(a)
    if iszero(a)
        T = eltype(a)
        if hasbasis(A)
            _coeff_elt_print(io, zero(T), first(basis(A)))
        else
            print(io, zero(T))
        end
    else #if hasbasis(A)
        #nzeros = findall(!iszero, coeffs(a))
        first = true
        for (idx, value) in zip(keys(coeffs(a)), values(coeffs(a)))
        #for (counter, idx) in enumerate(nzeros)
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
#    else
#        println(io, "algebra element without defined basis")
#        show(io, MIME("text/plain"), a.coeffs)
    end
end

function Base.show(io::IO, ::MIME"text/plain", mstr::LazyMStructure)
    print(io, "LazyMStructure of", object(basis(mstr)))
    if Base.haslength(basis(mstr))
        print(io, " and basis with $(length(basis(mstr))) elements")
    end
end
