# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module TestFusedArithmetic

using Test
using SparseArrays
import StarAlgebras as SA
import MutableArithmetics as MA

include(joinpath(@__DIR__, "..", "examples", "bivariate.jl"))
include(joinpath(@__DIR__, "..", "examples", "natural.jl"))
include(joinpath(@__DIR__, "..", "examples", "cheby.jl"))

SA.star(m::Monomial) = m

function sparse_coefficients(v)
    return SA.SparseCoefficients(findall(!iszero, v), filter(!iszero, v))
end

@testset "Fused algebra and term products" begin
    alg = SA.StarAlgebra(1.0, SA.FixedBasis(2.0 .^ (0:4)))
    for T in (Int, Rational{Int}, BigInt),
        storage in (identity, sparse, sparse_coefficients),
        op in (MA.add_mul, MA.sub_mul),
        left_term in (false, true),
        right_term in (false, true)

        f = SA.AlgebraElement(storage(T[1, -2, 1, 0, 0]), alg)
        a =
            left_term ? SA.Term(alg, 2, T(3)) :
            SA.AlgebraElement(storage(T[2, 3, 0, 0, 0]), alg)
        b =
            right_term ? SA.Term(alg, 2, T(2)) :
            SA.AlgebraElement(storage(T[-1, 1, 2, 0, 0]), alg)
        originals = MA.copy_if_mutable.((f, a, b))
        product =
            left_term ? (right_term ? T[0, 0, 6, 0, 0] : T[0, -3, 3, 6, 0]) :
            (right_term ? T[0, 4, 6, 0, 0] : T[-2, -1, 7, 6, 0])
        expected = SA.AlgebraElement(
            storage(MA.add_sub_op(op).(T[1, -2, 1, 0, 0], product)),
            alg,
        )
        buffer = MA.buffer_for(op, typeof(f), typeof(a), typeof(b))
        out = MA.mutable_copy(f)
        @test (@inferred MA.buffered_operate_to!!(buffer, out, op, f, a, b)) ===
              out
        @test out == expected
        @test (f, a, b) == originals
        @test (@inferred MA.buffered_operate!!(buffer, op, f, a, b)) === f
        @test f == expected
        @test (a, b) == originals[2:3]
        @test all(!iszero(last(p)) for p in SA.nonzero_pairs(SA.coeffs(f)))
        f = MA.mutable_copy(first(originals))
        @test MA.operate_to!(f, op, f, a, b) === f
        @test f == expected
    end

    for T in (Int, BigInt),
        storage in (identity, sparse, sparse_coefficients),
        op in (MA.add_mul, MA.sub_mul),
        alias in (:left, :right, :both)

        f = SA.AlgebraElement(storage(T[1, 2, 1, 0, 0]), alg)
        other = SA.AlgebraElement(storage(T[2, -1, 0, 0, 0]), alg)
        a = alias === :right ? other : f
        b = alias === :left ? other : f
        original = MA.mutable_copy(other)
        expected = op(f, a, b)
        @test (@inferred MA.operate!!(op, f, a, b)) === f
        @test f == expected
        @test other == original
    end

    f = SA.AlgebraElement(SA.SparseCoefficients([3], [6]), alg)
    a, b = SA.Term(alg, 2, 2), SA.Term(alg, 2, 3)
    @test iszero(MA.operate!!(MA.sub_mul, f, a, b))
    @test isempty(keys(SA.coeffs(f)))
    f = zero(Int, alg)
    r = @inferred MA.operate!!(MA.add_mul, f, a, SA.Term(alg, 2, 0.5))
    @test r !== f
    @test eltype(r) === Float64
    @test SA.coeffs(r)[3] == 1.0
    @test iszero(f)

    other_alg = SA.StarAlgebra(1.0, SA.FixedBasis(2.0 .^ (0:5)))
    for to_element in (identity, SA.algebra_element)
        a = to_element(SA.Term(alg, 2, 2))
        b = to_element(SA.Term(other_alg, 2, 3))
        @test_throws ArgumentError MA.operate!!(MA.add_mul, f, a, b)
        @test iszero(f)
    end
    a = SA.algebra_element(SA.Term(alg, 2, 2))
    b = SA.algebra_element(SA.Term(alg, 2, 3))
    originals = MA.mutable_copy.((a, b))
    @test_throws ArgumentError MA.operate_to!(a, MA.add_mul, f, a, b)
    @test_throws ArgumentError MA.operate_to!(b, MA.add_mul, f, a, b)
    @test (a, b) == originals
end

@testset "Fused scalar products" begin
    alg = SA.StarAlgebra(
        Monomial((0, 0)),
        SA.FixedBasis([Monomial((i, 0)) for i in 0:4]),
    )
    for T in (Int, Rational{Int}, BigInt),
        storage in (identity, sparse, sparse_coefficients),
        op in (MA.add_mul, MA.sub_mul),
        term in (false, true),
        left in (false, true)

        f = SA.AlgebraElement(storage(T[1, -2, 1, 0, 0]), alg)
        g =
            term ? SA.Term(alg, 2, T(3)) :
            SA.AlgebraElement(storage(T[0, 1, -3, 0, 2]), alg)
        # Numbers convert to T; an exactly representable Float64 must not widen it.
        a, b = left ? (2.0, g) : (g, 2.0)
        originals = MA.copy_if_mutable.((f, a, b))
        product = term ? T[0, 6, 0, 0, 0] : T[0, 2, -6, 0, 4]
        expected = SA.AlgebraElement(
            storage(MA.add_sub_op(op).(T[1, -2, 1, 0, 0], product)),
            alg,
        )
        if term
            @test a * b == SA.Term(alg, 2, T(6))
        end
        @test MA.promote_operation(op, typeof(f), typeof(a), typeof(b)) ===
              typeof(f)
        out = zero(f)
        buffer = MA.buffer_for(op, typeof(f), typeof(a), typeof(b))
        @test (@inferred MA.buffered_operate_to!!(buffer, out, op, f, a, b)) ===
              out
        @test out == expected
        @test (f, a, b) == originals
        @test (@inferred MA.buffered_operate!!(buffer, op, f, a, b)) === f
        @test f == expected
        @test (a, b) == originals[2:3]
        @test all(!iszero(last(p)) for p in SA.nonzero_pairs(SA.coeffs(f)))
        if T <: Integer
            @test_throws InexactError MA.operate!(op, f, 0.5, g)
            @test_throws InexactError MA.operate!(op, f, g, 0.5)
        end
    end

    for T in (Int, BigInt),
        storage in (identity, sparse, sparse_coefficients),
        op in (MA.add_mul, MA.sub_mul),
        left in (false, true)

        f = SA.AlgebraElement(storage(T[1, 2, 0, -1, 0]), alg)
        scalar = SA.coeffs(f)[2]
        a, b = left ? (scalar, f) : (f, scalar)
        expected = op(f, a, b)
        @test (@inferred MA.operate_to!!(f, op, f, a, b)) === f
        @test f == expected
        @test scalar == 2
    end

    other_alg = SA.StarAlgebra(
        Monomial((0, 0)),
        SA.FixedBasis([Monomial((i, 0)) for i in 0:5]),
    )
    for to_element in (identity, SA.algebra_element), left in (false, true)
        f = SA.algebra_element(SA.Term(alg, 1, 1))
        g = to_element(SA.Term(other_alg, 2, 3))
        a, b = left ? (2, g) : (g, 2)
        @test_throws ArgumentError MA.operate!!(MA.add_mul, f, a, b)
        @test_throws ArgumentError MA.operate_to!(f, MA.add_mul, f, a, b)
        @test SA.coeffs(f)[1] == 1
    end
    f = SA.algebra_element(SA.Term(alg, 1, 1))
    g = SA.algebra_element(SA.Term(alg, 2, 3))
    for (a, b) in ((2, g), (g, 2))
        @test_throws ArgumentError MA.operate_to!(g, MA.add_mul, f, a, b)
        @test SA.coeffs(g)[2] == 3
    end
    for g in (SA.Term(alg, 2, 0.5), SA.algebra_element(SA.Term(alg, 2, 0.5)))
        for (a, b) in ((2, g), (g, 2)), op in (MA.add_mul, MA.sub_mul)
            r = @inferred MA.operate!!(op, f, a, b)
            @test eltype(r) === Float64
            @test r == op(f, a, b)
            @test r !== f
            @test SA.coeffs(f)[1] == 1
        end
    end

    # Read-only or differently ordered coefficient inputs use canonicalization.
    for source in (
        SA.SparseCoefficients((3, 1, 1), (2, 3, -1)),
        SA.SparseCoefficients((3, 1), (2, 2), >),
    )
        f = SA.AlgebraElement(SA.SparseCoefficients([1, 2], [1, 4]), alg)
        g = SA.AlgebraElement(source, alg)
        @test MA.operate!(MA.sub_mul, f, 2, g) === f
        @test keys(SA.coeffs(f)) == [1, 2, 3]
        @test values(SA.coeffs(f)) == [-3, 4, -4]
    end
    f = SA.AlgebraElement(SA.SparseCoefficients([2, 1, 1], [3, 1, 1]), alg)
    @test MA.operate!(MA.add_mul, f, 2, f) === f
    @test keys(SA.coeffs(f)) == [1, 2]
    @test values(SA.coeffs(f)) == [6, 9]
end

@testset "Expanding products and coefficient order" begin
    alg = SA.StarAlgebra(ChebyPoly(0), ChebyMStruct(cheby_basis()))
    for to_element in (identity, SA.algebra_element),
        op in (MA.add_mul, MA.sub_mul)

        a = to_element(SA.Term(alg, 2, 3 // 1))
        b = to_element(SA.Term(alg, 3, 2 // 1))
        f = zero(Rational{Int}, alg)
        @test (@inferred MA.operate!!(op, f, a, b)) === f
        @test keys(SA.coeffs(f)) == [1, 5]
        @test values(SA.coeffs(f)) == fill(op === MA.add_mul ? 3 : -3, 2)
    end

    basis = SA.DiracBasis(["", "a", "b", "ab", "ba"])
    alg = SA.StarAlgebra("", SA.DiracMStructure(basis, *))
    A, B = [1 2; 0 1], [1 0; 3 1]
    for to_element in (identity, SA.algebra_element), reverse in (false, true)
        a = to_element(SA.Term(alg, "a", A))
        b = to_element(SA.Term(alg, "b", B))
        f = zero(Matrix{Int}, alg)
        left, right = reverse ? (b, a) : (a, b)
        @test MA.operate!!(MA.sub_mul, f, left, right) === f
        @test keys(SA.coeffs(f)) == [reverse ? "ba" : "ab"]
        @test only(values(SA.coeffs(f))) == -(reverse ? B * A : A * B)
        @test A == [1 2; 0 1]
        @test B == [1 0; 3 1]
    end
    for to_element in (identity, SA.algebra_element),
        op in (MA.add_mul, MA.sub_mul),
        left in (false, true)

        g = to_element(SA.Term(alg, "a", A))
        a, b = left ? (B, g) : (g, B)
        f = zero(Matrix{Int}, alg)
        @test (@inferred MA.operate!!(op, f, a, b)) === f
        @test only(values(SA.coeffs(f))) ==
              MA.add_sub_op(op)(left ? B * A : A * B)
        @test A == [1 2; 0 1]
        @test B == [1 0; 3 1]
    end
end

const biv_alg = bivariate_algebra()

function SA.term_product_style(
    ::typeof(SA.mstructure(biv_alg)),
    ::typeof(grlex),
)
    return SA.OrderedTermProduct()
end

function test_term_product_allocations(::Type{T}, op) where {T}
    a = SA.Term(biv_alg, (1, 0), T(2))
    b = SA.Term(biv_alg, (0, 1), T(3))
    f = zero(T, biv_alg)
    c = SA.coeffs(f)
    sizehint!(keys(c), 1)
    sizehint!(values(c), 1)
    @test (@inferred MA.operate!!(op, f, a, b)) === f
    @test keys(c) == [(1, 1)]
    @test values(c) == [op === MA.add_mul ? 6 : -6]
    MA.operate!(zero, f)
    @test (@allocated MA.operate!!(op, f, a, b)) == 0
    @test keys(c) == [(1, 1)]
    @test values(c) == [op === MA.add_mul ? 6 : -6]
    return
end

@testset "Allocation-free ordered term products" begin
    for T in (Int, Float64), op in (MA.add_mul, MA.sub_mul)
        test_term_product_allocations(T, op)
    end
end

function test_scalar_product_allocations(::Type{T}, op, left, term) where {T}
    t = SA.Term(biv_alg, (1, 0), T(2))
    g = term ? t : SA.algebra_element(t)
    a, b = left ? (T(3), g) : (g, T(3))
    f = zero(T, biv_alg)
    c = SA.coeffs(f)
    sizehint!(keys(c), 2)
    sizehint!(values(c), 2)
    @test (@inferred MA.operate!!(op, f, a, b)) === f
    MA.operate!(zero, f)
    @test (@allocated MA.operate!!(op, f, a, b)) == 0
    @test keys(c) == [(1, 0)]
    @test values(c) == [op === MA.add_mul ? 6 : -6]
    @test (@allocated MA.operate!!(op, f, a, b)) == 0
    @test values(c) == [op === MA.add_mul ? 12 : -12]
    return
end

@testset "Allocation-free scalar products" begin
    for T in (Int, Float64),
        op in (MA.add_mul, MA.sub_mul),
        left in (false, true),
        term in (false, true)

        test_scalar_product_allocations(T, op, left, term)
    end
end

end
