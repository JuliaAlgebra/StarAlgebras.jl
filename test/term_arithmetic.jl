# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

function SA.term_product_style(
    ::typeof(SA.mstructure(biv_alg)),
    ::typeof(grlex),
)
    return SA.OrderedTermProduct()
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

struct WeightedTermProduct{B} <:
       SA.MultiplicativeStructure{Monomial,NTuple{2,Int}}
    basis::B
end
function (ms::WeightedTermProduct)(
    a::NTuple{2,Int},
    b::NTuple{2,Int},
    ::Type{NTuple{2,Int}},
)
    return SA.SparseCoefficients((a .+ b,), (2,))
end
function SA.term_product_style(::WeightedTermProduct, ::typeof(grlex))
    return SA.OrderedTermProduct()
end

@testset "Product weights and coefficient order" begin
    alg = StarAlgebra(Monomial((0, 0)), WeightedTermProduct(biv_basis))
    A, B = [1 2; 0 1], [1 0; 3 1]
    original = deepcopy((A, B))
    t = Term(alg, (1, 0), A)
    g = SA.algebra_element(Term(alg, (0, 1), B))
    for left in (true, false)
        f = SA.algebra_element(Term(alg, (1, 1), zeros(Int, 2, 2)))
        a, b = left ? (t, g) : (g, t)
        @test MA.operate!!(MA.sub_mul, f, a, b) === f
        @test SA.coeffs(f)[(1, 1)] == -2 * (left ? A * B : B * A)
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
