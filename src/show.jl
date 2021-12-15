Base.show(io::IO, A::AbstractStarAlgebra) =
    print(io, "*-algebra of $(object(A))")
Base.show(io::IO, ::Type{<:StarAlgebra{O,T}}) where {O,T} =
    print(io, "StarAlgebra{$O, $T, …}")

__prints_with_minus(::Any) = false
__prints_with_minus(x::Real) = x < 0
__needs_parens(::Any) = false
__needs_parens(a::AlgebraElement) = true

function _coeff_elt_print(io, c, elt)
    print(io, c, "·")
    __needs_parens(elt) && print(io, "(")
    print(io, elt)
    __needs_parens(elt) && print(io, ")")
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
    elseif hasbasis(A)
        elts = String[]
        nzeros = findall(!iszero, coeffs(a))
        for (counter, idx) in enumerate(nzeros)
            c, elt = coeffs(a)[idx], basis(A)[idx]
            if counter == 1
                _coeff_elt_print(io, c, elt)
            else
                print(io, ' ')
                __prints_with_minus(c) || print(io, '+')
                _coeff_elt_print(io, c, elt)
                counter == length(nzeros) || print(io, ' ')
            end
        end
    else
        println(io, "algebra element without defined basis")
        show(io, MIME("text/plain"), a.coeffs)
    end
end

function Base.show(io::IO, ::MIME"text/plain", mstr::TrivialMStructure)
    Tw = _istwisted(mstr)
    l = length(basis(mstr))
    Tw && print(io, "twisted ")
    print(io, "TrivialMStructure over basis with $l elements")
    return
end
