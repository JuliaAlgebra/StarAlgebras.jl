abstract type MultiplicativeStructure{I} <: AbstractMatrix{I} end

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

struct TrivialMStructure{I,B<:AbstractBasis} <: MultiplicativeStructure{I}
    basis::B
end

TrivialMStructure(basis::AbstractBasis{T,I}) where {T,I} =
    TrivialMStructure{I,typeof(basis)}(basis)

basis(mstr::TrivialMStructure) = mstr.basis
Base.size(mstr::TrivialMStructure) = (l = length(basis(mstr)); (l, l))

Base.@propagate_inbounds function Base.getindex(
    mstr::TrivialMStructure,
    i::Integer,
    j::Integer,
)
    b = basis(mstr)
    g = i ≥ 0 ? b[i] : star(b[-i])
    h = j ≥ 0 ? b[j] : star(b[-j])
    gh = g * h
    gh in b || throw(ProductNotDefined(i, j, "$g · $h = $gh"))
    return basis(mstr)[gh]
end
