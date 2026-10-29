# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

using Test
using LinearAlgebra
import StarAlgebras as SA
import MutableArithmetics as MA

if !@isdefined bivariate_algebra
    include("../examples/bivariate.jl")
end

function _bivariate(alg, exponents, values)
    return SA.AlgebraElement(
        SA.SparseCoefficients(collect(exponents), collect(values), grlex),
        alg,
    )
end

@testset "Determinant" begin
    alg = bivariate_algebra()
    for T in (Bool, Int, BigInt)
        x = _bivariate(alg, [(1, 0)], T[1])
        y = _bivariate(alg, [(0, 1)], T[1])
        o = _bivariate(alg, [(0, 0)], T[1])
        z = _bivariate(alg, Tuple{Int,Int}[], T[])
        @test (@inferred det([x;;])) == x
        @test (@inferred det([x o; y o])) == x - y
        @test (@inferred det([x o; x o])) == z
        @test (@inferred det([x+y y; z y])) == (x + y) * y
        @test (@inferred det([x o o o; x y o o+o; z z z o; x z y z])) ==
              -x * y * y + (o + o) * x * y - x
        # The entries are not modified
        @test x == _bivariate(alg, [(1, 0)], T[1])
        @test o == _bivariate(alg, [(0, 0)], T[1])
    end
    for T in (Bool, Int, BigInt)
        x = SA.Term(alg, (1, 0), one(T))
        y = SA.Term(alg, (0, 1), one(T))
        o = SA.Term(alg, (0, 0), one(T))
        X = _bivariate(alg, [(1, 0)], T[1])
        Y = _bivariate(alg, [(0, 1)], T[1])
        O = _bivariate(alg, [(0, 0)], T[1])
        @test (@inferred det([x o; y o])) == X - Y
        @test (@inferred det([x;;])) == X
        @test (@inferred det([x o; y X + Y])) == X * X + X * Y - Y
        @test x == SA.Term(alg, (1, 0), one(T))
    end
    @test_throws DimensionMismatch det(fill(one(alg), 2, 3))
    @test_throws ArgumentError det(Matrix{typeof(one(alg))}(undef, 0, 0))
end
