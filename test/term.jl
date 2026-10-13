# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

using Test, LinearAlgebra, StarAlgebras
using SparseArrays
import StarAlgebras as SA
import MutableArithmetics as MA

if !@isdefined(bivariate_algebra)
    include("../examples/bivariate.jl")
end
if !@isdefined(NaturalNumbers)
    include("../examples/natural.jl")
end
if !@isdefined(cheby_basis)
    include("../examples/cheby.jl")
end

# Use the bivariate example algebra for testing
const biv_alg = bivariate_algebra()
const biv_basis = SA.basis(biv_alg)

# Helper: create a Term in the bivariate algebra
biv_term(coeff, exp::NTuple{2,Int}) = Term(biv_alg, exp, coeff)

function SA.term_product_style(
    ::typeof(SA.mstructure(biv_alg)),
    ::typeof(grlex),
)
    return SA.OrderedTermProduct()
end

SA.star(m::Monomial) = m
SA.star(t::Term{Int,typeof(biv_alg)}) = Term(t.algebra, t.index, coefficient(t))

function SA.promote_bases_with_maps(a::ChebyMStruct, b::ChebyMStruct)
    @assert SA.basis(a) === SA.basis(b)
    return (a, nothing), (b, nothing)
end

struct MatrixCoefficient
    value::Matrix{Int}
end
Base.zero(::Type{MatrixCoefficient}) = MatrixCoefficient(zeros(Int, 2, 2))
Base.iszero(x::MatrixCoefficient) = iszero(x.value)
function Base.:*(x::MatrixCoefficient, y::MatrixCoefficient)
    return MatrixCoefficient(x.value * y.value)
end

struct VectorTermStructure{B,S} <:
       SA.MultiplicativeStructure{Monomial,NTuple{2,Int}}
    basis::B
    style::S
end
function (ms::VectorTermStructure)(a, b, ::Type{NTuple{2,Int}})
    return SA.SparseCoefficients([a .+ b], [1])
end
SA.term_product_style(ms::VectorTermStructure, ::typeof(grlex)) = ms.style

function test_basis_weight_type(ms::M, ::Type{W}) where {M,W}
    @test (@inferred SA.value_type(M)) === W
    @test (@allocated SA.value_type(M)) == 0
    return
end

function test_product_promotion_allocations(a::A, b::B) where {A,B}
    MA.promote_operation(*, A, B)
    @test (@allocated MA.promote_operation(*, A, B)) == 0
    return
end

@testset "Term" begin
    @testset "Accessors" begin
        t = biv_term(3.0, (1, 2))
        @test coefficient(t) == 3.0
        @test basis_element(t) == Monomial((1, 2))
        @test parent(t) === biv_alg
    end

    @testset "iszero / zero" begin
        t = biv_term(3, (1, 0))
        @test !iszero(t)
        z = zero(t)
        @test iszero(z)
        @test coefficient(z) == 0
        @test basis_element(z) == Monomial((1, 0))
        @test iszero(SA.algebra_element(z))

        @test iszero(biv_term(0.0, (0, 0)))
        @test !iszero(biv_term(1, (0, 0)))
    end

    @testset "Coefficient mapping" begin
        alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0, 4.0]))
        for T in (Int, BigInt), nz in (false, true)
            for c in (
                T[1, 2, 3, 4],
                sparsevec(1:4, T[1, 2, 3, 4]),
                SA.SparseCoefficients(collect(1:4), T[1, 2, 3, 4]),
            )
                output = AlgebraElement(c, alg)
                t = Term(alg, 3, c[3])
                @test (@inferred SA.map_coefficients_to!(
                    output,
                    c -> 2c,
                    t;
                    nonzero = nz,
                )) === output
                @test [c[i] for i in 1:4] == [0, 0, 6, 0]
                @test coefficient(t) == 3
                @test SA.map_coefficients_to!(output, zero, t) === output
                @test iszero(output)
                @test isempty(collect(SA.nonzero_pairs(c)))
                @test SA.map_coefficients_to!(
                    output,
                    _ -> error("A zero term has no stored coefficients"),
                    zero(t);
                    nonzero = nz,
                ) === output
                @test iszero(output)
            end
        end
    end

    @testset "copy / mutable_copy" begin
        t = biv_term(5, (2, 3))
        t2 = copy(t)
        @test coefficient(t2) == 5
        @test basis_element(t2) == Monomial((2, 3))

        t3 = MA.mutable_copy(t)
        @test coefficient(t3) == coefficient(t)
        @test basis_element(t3) == basis_element(t)
    end

    @testset "copy constructor" begin
        t = biv_term(3, (1, 0))
        T = typeof(t)
        t2 = T(t)
        @test t2 === t
    end

    @testset "== / isequal" begin
        t1 = biv_term(3.0, (1, 0))
        t2 = biv_term(3.0, (1, 0))
        t3 = biv_term(4.0, (1, 0))
        t4 = biv_term(3.0, (0, 1))
        @test t1 == t2
        @test t1 != t3
        @test t1 != t4
        @test isequal(t1, t2)
        @test !isequal(t1, t3)
        # zero terms are equal regardless of index
        @test biv_term(0.0, (1, 0)) == biv_term(0.0, (0, 1))
        @test isequal(biv_term(0.0, (1, 0)), biv_term(0.0, (0, 1)))
    end

    @testset "hash" begin
        t1 = biv_term(3, (1, 0))
        t2 = biv_term(3, (1, 0))
        @test hash(t1) == hash(t2)
        # hash of zero term
        @test hash(biv_term(0, (1, 0))) == hash(0)
        # hash of term with coefficient 1 matches hash of basis_element
        t3 = biv_term(1, (2, 0))
        @test hash(t3) == hash(Monomial((2, 0)))
    end

    @testset "convert" begin
        t = biv_term(3, (1, 0))
        A = typeof(biv_alg)
        I = NTuple{2,Int}
        t2 = convert(Term{Float64,A,I}, t)
        @test coefficient(t2) == 3.0
        @test basis_element(t2) == Monomial((1, 0))
        # identity convert
        t3 = convert(typeof(t), t)
        @test t3 === t
        @test convert(Term{Number,A,I}, t) isa Term{Number,A,I}

        ae = AlgebraElement(
            SA.SparseCoefficients((t.index,), (coefficient(t),), grlex),
            biv_alg,
        )
        @test convert(typeof(ae), ae) === ae
        for T in (Rational{Int}, Float64)
            R = SA.similar_type(typeof(SA.algebra_element(t)), T)
            converted = convert(R, ae)
            @test converted isa R
            @test parent(converted) === biv_alg
            @test SA.coeffs(converted)[t.index] == T(3)
            @test SA.coeffs(converted).isless === grlex
            @test convert(R, converted) === converted
            @test convert(R, t) isa R
            @test convert(R, t) == converted
            @test iszero(convert(R, zero(t)))
        end
        fractional = SA.algebra_element(biv_term(1 // 2, t.index))
        @test_throws InexactError convert(
            typeof(SA.algebra_element(t)),
            fractional,
        )
    end

    @testset "broadcastable / ndims" begin
        t = biv_term(2, (1, 0))
        @test ndims(t) == 0
        @test ndims(typeof(t)) == 0
        @test Base.broadcastable(t) isa Ref
    end

    @testset "negation" begin
        t = biv_term(3, (1, 0))
        @test -t == biv_term(-3, (1, 0))
    end

    @testset "Scalar multiplication and promotion" begin
        t = biv_term(2, (1, 0))
        u = biv_term(3.0, (0, 1))
        for x in (t, SA.algebra_element(t))
            @test +x === x
            @test *(x) === x
        end
        for (a, expected) in
            ((2.0, 4.0), (1.5, 3.0), (1 // 2, 1 // 1), (1 + im, 2 + 2im))
            for args in ((a, t), (t, a))
                result = @inferred *(args...)
                @test result == biv_term(expected, (1, 0))
                @test coefficient(result) === expected
                @test (@inferred MA.promote_operation(*, typeof.(args)...)) ===
                      typeof(result)
            end
        end
        @test (@inferred 1.5 * u) == biv_term(4.5, (0, 1))
        @test (@inferred u * 1.5) == biv_term(4.5, (0, 1))
        @test_throws MethodError "x" * t
        @test_throws MethodError t * "x"

        a = MatrixCoefficient([1 2; 0 1])
        t_matrix = biv_term(MatrixCoefficient([1 0; 3 1]), (1, 0))
        @test coefficient(a * t_matrix).value == [7 2; 3 1]
        @test coefficient(t_matrix * a).value == [1 2; 3 7]

        p, q = SA.algebra_element(t), SA.algebra_element(u)
        for (a, b, T) in
            ((t, u, typeof(u)), (p, q, typeof(q)), (p, u, typeof(q)))
            @test (@inferred Base.vect(a, b)) isa Vector{T}
            @test (@inferred Base.vect(b, a)) isa Vector{T}
            @test [a, b] == [convert(T, a), convert(T, b)]
        end
    end

    @testset "Term + Term → AlgebraElement" begin
        t1 = biv_term(3, (2, 0))
        t2 = biv_term(2, (0, 3))
        ae = t1 + t2
        @test ae isa SA.AlgebraElement
        # Check it has two terms
        c = SA.coeffs(ae)
        @test length(collect(SA.keys(c))) == 2
        @test c.isless === grlex
        @test collect(keys(c)) == [(2, 0), (0, 3)]

        # Same index: should combine
        t3 = biv_term(5, (2, 0))
        ae2 = t1 + t3
        c2 = SA.coeffs(ae2)
        @test length(collect(SA.keys(c2))) == 1
    end

    @testset "Term - Term → AlgebraElement" begin
        t1 = biv_term(3, (1, 0))
        t2 = biv_term(2, (1, 0))
        ae = t1 - t2
        @test ae isa SA.AlgebraElement
        c = SA.coeffs(ae)
        @test length(collect(SA.keys(c))) == 1
        @test iszero(t1 - t1)
    end

    @testset "Ordered term products stay terms" begin
        test_basis_weight_type(SA.mstructure(biv_alg), Int)
        for T in (Int, Float64, Rational{Int}, BigInt)
            a, b = biv_term(T(2), (1, 0)), biv_term(T(3), (0, 2))
            product = @inferred a * b
            @test product isa Term{T}
            @test (@inferred MA.promote_operation(*, typeof(a), typeof(b))) ===
                  typeof(product)
            @test product == biv_term(T(6), (1, 2))
            z = @inferred zero(a) * b
            @test z isa typeof(product)
            @test iszero(z)
            @test a == biv_term(T(2), (1, 0))
            @test b == biv_term(T(3), (0, 2))
            if T === Int
                test_product_promotion_allocations(a, b)
                @test (@allocated a * b) == 0
            end
        end
        # The trait decides the representation even for vector-backed products.
        for (style, R) in (
            (SA.OrderedTermProduct(), Term{Int}),
            (SA.GeneralTermProduct(), AlgebraElement{Int}),
        )
            alg = StarAlgebra(
                Monomial((0, 0)),
                VectorTermStructure(biv_basis, style),
            )
            test_basis_weight_type(SA.mstructure(alg), Int)
            a, b = Term(alg, (1, 0), 3), Term(alg, (0, 1), 1)
            product = @inferred a * b
            @test product isa R
            expected = Term(alg, (1, 1), 3)
            @test product == (
                style isa SA.OrderedTermProduct ? expected :
                SA.algebra_element(expected)
            )
            @test (@inferred MA.promote_operation(*, typeof(a), typeof(b))) ===
                  typeof(product)
            @test (@inferred zero(a) * b) isa typeof(product)
        end

        alg = StarAlgebra("", SA.DiracBasis(["", "a", "b", "ab", "ba"]))
        a, b = Term(alg, "a", 2), Term(alg, "b", 3)
        @test a * b isa AlgebraElement{Int}
        @test (@inferred a * b) == SA.algebra_element(Term(alg, "ab", 6))
        @test (@inferred b * a) == SA.algebra_element(Term(alg, "ba", 6))
    end

    @testset "Term product can expand" begin
        alg = SA.StarAlgebra(ChebyPoly(0), ChebyMStruct(cheby_basis()))
        a, b = Term(alg, 2, 3 // 1), Term(alg, 3, 2 // 1)
        product = @inferred a * b
        @test product isa AlgebraElement{Rational{Int}}
        @test collect(keys(SA.coeffs(product))) == [1, 5]
        @test collect(values(SA.coeffs(product))) == [3 // 1, 3 // 1]
        @test iszero(@inferred zero(a) * b)
        @test SA.coeffs(@inferred a * a)[0] == 9 // 2
    end

    @testset "Basis weight promotion" begin
        ms = ChebyMStruct(cheby_basis())
        for structure in (ms, SA.MTable(ms, (0, 0)))
            test_basis_weight_type(structure, Rational{Int})
            alg = StarAlgebra(ChebyPoly(0), structure)
            for (T, R) in (
                (Int, Rational{Int}),
                (BigInt, Rational{BigInt}),
                (Float64, Float64),
            )
                t = Term(alg, 2, T(3))
                p = SA.algebra_element(t)
                for (a, b) in ((t, t), (t, p), (p, t), (p, p))
                    product = @inferred a * b
                    @test product isa AlgebraElement{R}
                    @test (@inferred MA.promote_operation(
                        *,
                        typeof(a),
                        typeof(b),
                    )) === typeof(product)
                    @test collect(keys(SA.coeffs(product))) == [0, 4]
                    @test collect(values(SA.coeffs(product))) == R[9//2, 9//2]
                    z = @inferred zero(a) * b
                    @test z isa typeof(product)
                    @test iszero(z)
                    if T === Int
                        test_product_promotion_allocations(a, b)
                    end
                end
                @test coefficient(t) == 3
                @test SA.coeffs(p)[2] == 3
                # Scalar multiplication does not multiply basis elements.
                @test (@inferred 2 * p) isa AlgebraElement{T}
                @test (@inferred p * 2) isa AlgebraElement{T}
            end
        end
    end

    @testset "AlgebraElement coefficient scaling" begin
        p = SA.algebra_element(biv_term(2.5, (1, 0)))
        for scaled in (2 * p, p * 2)
            @test eltype(scaled) === Float64
            @test SA.coeffs(scaled)[(1, 0)] == 5.0
        end
        @test iszero(0 * p)
        @test_throws MethodError "x" * p
        @test_throws MethodError p * "x"
        p_int = SA.algebra_element(biv_term(2, (1, 0)))
        for input in (p_int, SA.coeffs(p_int))
            for (op, args, expected) in (
                (*, (0.5, input), 1.0),
                (*, (input, 0.5), 1.0),
                (/, (input, 1.5), 2 / 1.5),
                (//, (input, 3 // 2), 4 // 3),
                (div, (input, 1.5), 1.0),
            )
                result = @inferred op(args...)
                c = result isa AlgebraElement ? SA.coeffs(result) : result
                @test c[(1, 0)] === expected
                if result isa AlgebraElement
                    @test (@inferred MA.promote_operation(
                        op,
                        typeof.(args)...,
                    )) === typeof(result)
                end
            end
        end
        for args in ((0.5, p_int), (p_int, 0.5))
            output = similar(p_int, Float64)
            @test (@inferred MA.operate_to!(output, *, args...)) === output
            @test output == SA.algebra_element(biv_term(1.0, (1, 0)))
        end
        @test_throws InexactError MA.operate_to!(similar(p_int), *, p_int, 0.75)
        @test SA.coeffs(p_int)[(1, 0)] === 2

        a = MatrixCoefficient([1 2; 0 1])
        p_matrix =
            SA.algebra_element(biv_term(MatrixCoefficient([1 0; 3 1]), (1, 0)))
        @test SA.coeffs(a * p_matrix)[(1, 0)].value == [7 2; 3 1]
        @test SA.coeffs(p_matrix * a)[(1, 0)].value == [1 2; 3 7]
    end

    @testset "Term + AlgebraElement" begin
        t = biv_term(3, (1, 0))
        ae = t + t
        @test ae isa SA.AlgebraElement
        t2 = biv_term(2, (0, 1))
        ae2 = t2 + ae
        @test ae2 isa SA.AlgebraElement
        ae3 = ae + t2
        @test ae3 isa SA.AlgebraElement
    end

    @testset "algebra_element(::Term)" begin
        t = biv_term(3, (1, 0))
        ae = SA.algebra_element(t)
        @test ae isa SA.AlgebraElement
        @test parent(ae) === biv_alg
    end

    @testset "various coefficient types" begin
        # Float64
        t = biv_term(2.5, (1, 1))
        @test coefficient(t) == 2.5
        @test !iszero(t)

        # Complex
        tc = biv_term(1 + 2im, (0, 1))
        @test coefficient(tc) == 1 + 2im
        @test !iszero(tc)
        @test iszero(zero(tc))
        @test coefficient(zero(tc)) == 0 + 0im
        @test LinearAlgebra.dot(tc, tc) == 5

        # Rational
        tr = biv_term(3 // 4, (2, 0))
        @test coefficient(tr) == 3 // 4
        @test coefficient(zero(tr)) == 0 // 1
    end
end
