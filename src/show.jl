function Base.show(io::IO, A::AbstractStarAlgebra)
    ioc = IOContext(io, :limit => true, :compact => true)
    return print(ioc, "*-algebra of ", object(A))
end

__needs_parens(::Any) = false
__needs_parens(a::AlgebraElement) = true

# `print_coefficient` is inspired from MultivariatePolynomials.jl

# `Int`, `Float64` don't support MIME"text/latex".
# We could add a check with `showable` if a `Real` subtype supports it and
# the feature is requested.
print_coefficient(io::IO, ::MIME, coeff::Real) = print(io, coeff)
# Scientific notation does not display well in LaTeX so we rewrite it
function print_coefficient(io::IO, ::MIME"text/latex", coeff::AbstractFloat)
    s = string(coeff)
    if occursin('e', s)
        s = replace(s, 'e' => " \\cdot 10^{") * '}'
    end
    return print(io, s)
end

trim_LaTeX(_, s::AbstractString) = s

function trim_LaTeX(::MIME"text/latex", s::AbstractString)
    i = firstindex(s)
    j = lastindex(s)
    while true
        if i < j && isspace(s[i])
            i = nextind(s, i)
        elseif i < j && isspace(s[j])
            j = prevind(s, j)
        elseif i < j && s[i] == '$' && s[j] == '$'
            i = nextind(s, i)
            j = prevind(s, j)
        elseif i < j && (
            (s[i:nextind(s, i)] == "\\(" && s[prevind(s, j):j] == "\\)") ||
            (s[i:nextind(s, i)] == "\\[" && s[prevind(s, j):j] == "\\]")
        )
            i = nextind(s, i, 2)
            j = prevind(s, j, 2)
        else
            return s[i:j]
        end
    end
end

# JuMP expressions supports LaTeX output so `showable` will return `true`
# for them. It is important for anonymous variables to display properly as well:
# https://github.com/jump-dev/SumOfSquares.jl/issues/256
# Since they add `$$` around it, we need to trim it with `trim_LaTeX`
function print_coefficient(io::IO, mime, coeff)
    print(io, "(")
    print_mime(io, mime, coeff)
    print(io, ")")
    return
end

function print_mime(io::IO, mime, x)
    if showable(mime, x)
        print(io, trim_LaTeX(mime, sprint(show, mime, x)))
    else
        show(io, x)
    end
end

isnegative(x::Real) = x < 0
isnegative(x) = false

_print_dot(io, ::MIME"text/latex") = print(io, " \\cdot ")
_print_dot(io, ::MIME) = print(io, '·')

function _coeff_elt_print(io, mime, c, elt)
    print_coefficient(io, mime, c)
    _print_dot(io, mime)
    __needs_parens(elt) && print(io, '(')
    print_mime(io, mime, elt)
    __needs_parens(elt) && print(io, ')')
    return
end

Base.print(io::IO, a::AlgebraElement) = show(io, MIME"text/print"(), a)
Base.show(io::IO, a::AlgebraElement) = show(io, MIME"text/plain"(), a)

function Base.show(io::IO, mime::MIME"text/latex", a::AlgebraElement)
    print(io, "\$\$ ")
    _show(io, mime, a)
    return print(io, " \$\$")
end

# If the MIME is not specified, IJulia thinks that it supports images, ...
# and then use the result of show and tries to interpret it as an svg, ...
# We need the two methods to avoid ambiguity
function Base.show(io::IO, mime::MIME"text/plain", a::AlgebraElement)
    return _show(io, mime, a)
end
function Base.show(io::IO, mime::MIME"text/print", a::AlgebraElement)
    return _show(io, mime, a)
end

function _show(io::IO, mime, a::AlgebraElement)
    A = parent(a)
    if iszero(a)
        T = value_type(coeffs(a))
        _coeff_elt_print(io, mime, zero(T), first(basis(A)))
    else
        _first = true
        for (idx, value) in nonzero_pairs(coeffs(a))
            c, elt = value, basis(A)[idx]
            if _first
                _coeff_elt_print(io, mime, c, elt)
                _first = false
            else
                neg = isnegative(c)
                print(io, ' ', neg ? '-' : '+', ' ')
                _coeff_elt_print(io, mime, neg ? abs(c) : c, elt)
            end
        end
    end
end
