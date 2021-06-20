abstract type MultiplicativeStructure{Twisted,I} <: AbstractMatrix{I} end

_istwisted(::MultiplicativeStructure{T}) where {T} = T
_product(ms::MultiplicativeStructure, g, h) = _product(Val(_istwisted(ms)), g, h)

_product(::Val{false}, g, h) = g * h
_product(::Val{true}, g, h) = star(g) * h

struct ProductNotDefined <: Exception
    i::Any
    j::Any
    msg::Any
end

function Base.showerror(io::IO, ex::ProductNotDefined)
    print(io, "Product of elements $(ex.i) and $(ex.j) is not defined on the basis")
    print(io, " or the multiplicative structure could not be completed")
    if isdefined(ex, :msg)
        print(io, ": $(ex.msg)")
    end
    print(io, ".")
end

struct TrivialMStructure{Tw,I,B<:AbstractBasis} <: MultiplicativeStructure{Tw,I}
    basis::B
end

TrivialMStructure{Tw}(basis::AbstractBasis{T,I}) where {Tw,T,I} =
    TrivialMStructure{Tw,I,typeof(basis)}(basis)

basis(mstr::TrivialMStructure) = mstr.basis
Base.size(mstr::TrivialMStructure) = (l = length(basis(mstr)); (l, l))

Base.@propagate_inbounds function Base.getindex(mstr::TrivialMStructure, i::Integer, j::Integer)
    @boundscheck checkbounds(mstr, i, j)
    b = basis(mstr)
    g, h = b[i], b[j]
    gh = _product(mstr, g, h)
    gh in b || throw(ProductNotDefined(i, j, "$g · $h = $gh"))
    return b[gh]
end    
