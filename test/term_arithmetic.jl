# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

using Test, StarAlgebras
import StarAlgebras as SA
import MutableArithmetics as MA

if !@isdefined(biv_alg)
    include("term.jl")
end

@testset "Term products in a finite basis" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 4.0]))
    t = Term(alg, 2, 2)
    for c in (
        [1, 3, 0],
        sparsevec([1, 2], [1, 3], 3),
        SA.SparseCoefficients([1, 2], [1, 3]),
    )
        p = AlgebraElement(c, alg)
        for (a, b) in ((t, p), (p, t))
            product = @inferred a * b
            @test [SA.coeffs(product)[i] for i in 1:3] == [0, 2, 6]
            @test typeof(product) ===
                  MA.promote_operation(*, typeof(a), typeof(b))
            @test [SA.coeffs(p)[i] for i in 1:3] == [1, 3, 0]
            @test coefficient(t) == 2
        end
    end
end

mutable struct CountingTermOrder
    comparisons::Int
end
function (lt::CountingTermOrder)(a, b)
    lt.comparisons += 1
    return grlex(a, b)
end
function SA.term_product_style(
    ::typeof(SA.mstructure(biv_alg)),
    ::CountingTermOrder,
)
    return SA.OrderedTermProduct()
end

@testset "Ordered term accumulation" begin
    ms = SA.mstructure(biv_alg)
    @test SA.term_product_style(ms, grlex) isa SA.OrderedTermProduct
    @test SA.term_product_style(ms, isless) isa SA.GeneralTermProduct
    for T in (Int, Rational{Int}, BigInt),
        op in (MA.add_mul, MA.sub_mul),
        left in (true, false)

        t = biv_term(T(2), (1, 0))
        g = biv_term(T(3), (0, 0)) + biv_term(T(1), (2, 0))
        f = biv_term(T(6), (1, 0)) + biv_term(T(1), (0, 2))
        original = deepcopy((t, g))
        a, b = left ? (t, g) : (g, t)
        expected = op(f, a, b)
        buffer = MA.buffer_for(op, typeof(f), typeof(a), typeof(b))
        @test buffer === nothing
        out = zero(f)
        f_original = deepcopy(f)
        @test MA.buffered_operate_to!!(buffer, out, op, f, a, b) === out
        @test out == expected
        @test f == f_original
        @test MA.buffered_operate!!(buffer, op, f, a, b) === f
        @test f == expected
        @test (t, g) == original
        @test all(!iszero, values(SA.coeffs(f)))
        @test MA.operate!!(op, f, zero(t), g) === f
        @test f == expected
        @test MA.operate!!(op, f, t, zero(g)) === f
        @test f == expected

        # Backwards merging must preserve unread entries even when f === g.
        expected = left ? op(f, t, f) : op(f, f, t)
        a, b = left ? (t, f) : (f, t)
        @test MA.buffered_operate!(buffer, op, f, a, b) === f
        @test f == expected
    end

    lt = CountingTermOrder(0)
    c = SA.SparseCoefficients([(i, 0) for i in 0:99], ones(Int, 100), lt)
    f, g = AlgebraElement(copy(c), biv_alg), AlgebraElement(copy(c), biv_alg)
    t = biv_term(1, (1, 0))
    buffer = MA.buffer_for(MA.sub_mul, typeof(f), typeof(t), typeof(g))
    @test buffer === nothing
    @test MA.buffered_operate!(buffer, MA.sub_mul, f, t, g) === f
    @test lt.comparisons <= 200
    @test collect(SA.nonzero_pairs(SA.coeffs(f))) ==
          [(0, 0) => 1, (100, 0) => -1]
    @test MA.buffered_operate!(buffer, MA.add_mul, f, t, g) === f
    @test keys(SA.coeffs(f)) == keys(c)
    @test values(SA.coeffs(f)) == values(c)

    # Different storage orders must use the general accumulation path.
    f = AlgebraElement(SA.SparseCoefficients([(0, 2)], [1], grlex), biv_alg)
    g = AlgebraElement(SA.SparseCoefficients([(0, 2), (1, 0)], [1, 2]), biv_alg)
    expected = MA.sub_mul(f, t, g)
    @test MA.operate!!(MA.sub_mul, f, t, g) === f
    @test f == expected
end

@testset "Accumulation of two terms" begin
    a, b = biv_term(2, (1, 0)), biv_term(3, (0, 1))
    for order in (grlex, isless), op in (MA.add_mul, MA.sub_mul)
        f = AlgebraElement(SA.SparseCoefficients([(0, 0)], [1], order), biv_alg)
        original = copy(f)
        expected = op(f, a, b)
        output = zero(f)
        @test (@inferred MA.operate_to!!(output, op, f, a, b)) === output
        @test output == expected
        @test f == original
        @test (@inferred MA.operate!!(op, f, a, b)) === f
        @test f == expected
        @test a == biv_term(2, (1, 0))
        @test b == biv_term(3, (0, 1))
    end
end

struct WeightedTermProduct{B,W} <:
       SA.MultiplicativeStructure{Monomial,NTuple{2,Int}}
    basis::B
    weight::W
end
SA.value_type(::Type{<:WeightedTermProduct{B,W}}) where {B,W} = W
function (ms::WeightedTermProduct)(
    a::NTuple{2,Int},
    b::NTuple{2,Int},
    ::Type{NTuple{2,Int}},
)
    return SA.SparseCoefficients((a .+ b,), (ms.weight,))
end
function SA.term_product_style(::WeightedTermProduct, ::typeof(grlex))
    return SA.OrderedTermProduct()
end

@testset "Product weights and coefficient order" for weight in (2, 1 // 2)
    alg = StarAlgebra(Monomial((0, 0)), WeightedTermProduct(biv_basis, weight))
    A, B = [1 2; 0 1], [1 0; 3 1]
    original = deepcopy((A, B))
    t = Term(alg, (1, 0), A)
    g = SA.algebra_element(Term(alg, (0, 1), B))
    u = Term(alg, (0, 1), B)
    for (a, b, expected) in ((t, u, weight * (A * B)), (u, t, weight * (B * A)))
        product = @inferred a * b
        @test product isa Term{typeof(expected)}
        @test product.index == (1, 1)
        @test coefficient(product) == expected
        @test (@inferred MA.promote_operation(*, typeof(a), typeof(b))) ===
              typeof(product)
    end
    for left in (true, false)
        expected = -weight * (left ? A * B : B * A)
        f = SA.algebra_element(Term(alg, (1, 1), zero(expected)))
        a, b = left ? (t, g) : (g, t)
        product = @inferred a * b
        @test (@inferred MA.promote_operation(*, typeof(a), typeof(b))) ===
              typeof(product)
        @test SA.coeffs(product)[(1, 1)] == -expected
        @test (@inferred MA.operate!!(MA.sub_mul, f, a, b)) === f
        @test SA.coeffs(f)[(1, 1)] == expected
        @test (coefficient(t), SA.coeffs(g)[(0, 1)]) == original
    end
end

@testset "Expanding term products" begin
    alg = StarAlgebra(ChebyPoly(0), ChebyMStruct(cheby_basis()))
    @test SA.term_product_style(SA.mstructure(alg), isless) isa
          SA.GeneralTermProduct
    t = Term(alg, 2, 3 // 1)
    for op in (MA.add_mul, MA.sub_mul), left in (true, false)
        g = SA.algebra_element(Term(alg, 3, 2 // 1))
        expected = left ? op(g, t, g) : op(g, g, t)
        a, b = left ? (t, g) : (g, t)
        @test MA.operate!!(op, g, a, b) === g
        @test g == expected
    end
    t = Term(alg, 2, 3)
    g = SA.algebra_element(t)
    for op in (MA.add_mul, MA.sub_mul), (a, b) in ((t, g), (g, t), (g, g))
        expected = op(g, a, b)
        @test (@inferred MA.promote_operation(
            op,
            typeof(g),
            typeof(a),
            typeof(b),
        )) === typeof(expected)
        original = copy(g)
        output = @inferred MA.operate!!(op, g, a, b)
        @test output isa AlgebraElement{Rational{Int}}
        @test output == expected
        @test g == original
        output = zero(expected)
        @test (@inferred MA.operate_to!(output, op, g, a, b)) === output
        @test output == expected
    end
end

word_order(a, b) = isless((length(a), a), (length(b), b))
function SA.term_product_style(
    ::SA.DiracMStructure{String,String,B,typeof(*)},
    ::typeof(word_order),
) where {B}
    return SA.OrderedTermProduct()
end

@testset "Noncommutative ordered products" begin
    basis = SA.DiracBasis(["", "a", "b", "ab", "ba", "aba", "baa"])
    alg = StarAlgebra("", SA.DiracMStructure(basis, *))
    t = Term(alg, "a", 2)
    g = AlgebraElement(
        SA.SparseCoefficients(["", "b", "ba"], [1, 2, 3], word_order),
        alg,
    )
    for left in (true, false)
        f = AlgebraElement(
            SA.SparseCoefficients(String[], Int[], word_order),
            alg,
        )
        a, b = left ? (t, g) : (g, t)
        @test MA.operate!!(MA.sub_mul, f, a, b) === f
        @test keys(SA.coeffs(f)) ==
              (left ? ["a", "ab", "aba"] : ["a", "ba", "baa"])
        @test values(SA.coeffs(f)) == [-2, -4, -6]
    end
end
