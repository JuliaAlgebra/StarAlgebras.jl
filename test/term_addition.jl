# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module TestTermAddition

using Test
using SparseArrays
import StarAlgebras as SA
import MutableArithmetics as MA

@testset "In-place term addition" begin
    alg = SA.StarAlgebra(1, SA.FixedBasis([1, 2, 3, 4]))
    for T in (Int, Float64, BigInt)
        for c in (
            T[1, 0, 3, 0],
            sparsevec([1, 3], T[1, 3], 4),
            SA.SparseCoefficients([1, 3], T[1, 3]),
        )
            f = SA.AlgebraElement(c, alg)
            for (index, value, expected) in (
                (3, -3, [1, 0, 0, 0]),
                (2, 2, [1, 2, 0, 0]),
                (4, 4, [1, 2, 0, 4]),
                (1, 0, [1, 2, 0, 4]),
            )
                t = SA.Term(alg, index, T(value))
                @test (@inferred MA.add!!(f, t)) === f
                @test [c[i] for i in 1:4] == expected
                @test SA.coefficient(t) == value
            end
            if c isa Union{SparseVector,SA.SparseCoefficients}
                @test [i => v for (i, v) in SA.nonzero_pairs(c)] == [1 => 1, 2 => 2, 4 => 4]
            end
        end
    end

    f = SA.AlgebraElement(sparsevec([1], [1], 4), alg)
    t = SA.Term(alg, 2, 0.5)
    r = @inferred MA.add!!(f, t)
    @test r !== f
    @test eltype(r) === Float64
    @test SA.coeffs(r) == [1.0, 0.5, 0, 0]
    @test SA.coeffs(f) == [1, 0, 0, 0]

    other_alg = SA.StarAlgebra(1, SA.FixedBasis([1, 2, 3, 4, 5]))
    t = SA.Term(other_alg, 5, 2)
    @test MA.mutability(f, +, f, t) isa MA.IsMutable
    @test_throws ArgumentError MA.add!!(f, t)
    @test SA.coeffs(f) == [1, 0, 0, 0]

    for c in (
        BigInt[0, 0, 0, 0],
        spzeros(BigInt, 4),
        SA.SparseCoefficients(Int[], BigInt[]),
    )
        f = SA.AlgebraElement(c, alg)
        t = SA.Term(alg, 2, big(3))
        @test MA.add!!(f, t) === f
        MA.operate!(+, c[2], 1)
        @test c[2] == 4
        @test SA.coefficient(t) == 3
    end
end

function test_allocations(::Type{T}, n, cancel) where {T}
    alg = SA.StarAlgebra(1, SA.FixedBasis(collect(1:(n+1))))
    c = SA.SparseCoefficients(collect(1:n), ones(T, n))
    f = SA.AlgebraElement(c, alg)
    t = SA.Term(alg, cancel && n > 0 ? 1 : n + 1, cancel ? -one(T) : one(T))
    expected = [ones(T, n); zero(T)]
    expected[t.index] += SA.coefficient(t)
    sizehint!(keys(c), n + 1)
    sizehint!(values(c), n + 1)
    @test (@inferred MA.add!!(f, t)) === f
    @test [c[i] for i in 1:(n+1)] == expected
    resize!(keys(c), n)
    resize!(values(c), n)
    copyto!(keys(c), 1:n)
    fill!(values(c), one(T))
    @test (@allocated MA.add!!(f, t)) == 0
    @test [c[i] for i in 1:(n+1)] == expected
    return
end

@testset "Allocation-free term addition" begin
    for T in (Int, Float64), n in (0, 1, 100), cancel in (false, true)
        test_allocations(T, n, cancel)
    end
end

end
