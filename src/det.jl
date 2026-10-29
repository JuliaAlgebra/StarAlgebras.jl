# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

"""
    LinearAlgebra.det(M::AbstractMatrix{<:Union{Term,AlgebraElement}})

Determinant of the square matrix `M` computed by cofactor expansion along the
first column. The entries of `M` should commute with each other as the
determinant is not well defined otherwise.
"""
function LinearAlgebra.det(M::AbstractMatrix{<:Union{Term,AlgebraElement}})
    n = LinearAlgebra.checksquare(M)
    if iszero(n)
        msg = "the algebra of an empty matrix is unknown"
        throw(ArgumentError("Cannot compute the determinant: $msg."))
    end
    T = eltype(M)
    R = MA.promote_operation(+, MA.promote_operation(*, T, T), T)
    return _cofactor_det(R, M, collect(axes(M, 1)), firstindex(M, 2))
end

# Determinant of `M[rows, col:end]`
function _cofactor_det(::Type{R}, M, rows, col) where {R}
    total = MA.zero_for(R, (M[first(rows), col],))
    if isone(length(rows))
        return MA.operate!!(+, total, M[only(rows), col])
    end
    for (k, i) in enumerate(rows)
        a = M[i, col]
        iszero(a) && continue
        minor = _cofactor_det(R, M, deleteat!(copy(rows), k), col + 1)
        op = isodd(k) ? MA.add_mul : MA.sub_mul
        total = MA.operate!!(op, total, a, minor)
    end
    return total
end
