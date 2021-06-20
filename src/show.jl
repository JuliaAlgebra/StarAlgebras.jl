Base.show(io::IO, A::AbstractStarAlgebra) = print(io, "*-Algebra of $(object(A))")

__prints_with_minus(x) = false
__prints_with_minus(x::Real) = x < 0

function Base.show(io::IO, a::AlgebraElement)
    A = parent(a)
    if iszero(a)
        T = eltype(a)
        print(io, "$(zero(T))*$(one(object(A)))")
    elseif hasbasis(A)
        elts = String[]
        nzeros = findall(!iszero, coeffs(a))
        for (counter, idx) in enumerate(nzeros)
            c, elt = coeffs(a)[idx], basis(A)[idx]
            if counter == 1
                print(io, c, '·', elt)
                length(nzeros) > 1 && print(io, ' ')
            else
                __prints_with_minus(c) || print(io, '+')
                print(io, c, '·', elt)
                counter == length(nzeros) || print(io, ' ')
            end
        end
    else
        println(io, "Algebra element without defined basis")
        show(io, MIME("text/plain"), a.coeffs)
    end
end

function Base.show(io::IO, ::MIME"text/plain", mstr::TrivialMStructure)
    Tw = _istwisted(mstr)
    l = length(basis(mstr))
    print(io, "TrivialMStructure{$Tw} over basis with $l elements")
end
