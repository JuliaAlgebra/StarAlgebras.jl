# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# Use the bivariate example algebra for testing
const biv_alg = bivariate_algebra()
const biv_basis = SA.basis(biv_alg)

# Helper: create a Term in the bivariate algebra
biv_term(coeff, exp::NTuple{2,Int}) = Term(biv_alg, exp, coeff)

SA.star(m::Monomial) = m
SA.star(t::Term{Int,typeof(biv_alg)}) = Term(t.algebra, t.index, coefficient(t))

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

        @test iszero(biv_term(0.0, (0, 0)))
        @test !iszero(biv_term(1, (0, 0)))
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

    @testset "Term + Term → AlgebraElement" begin
        t1 = biv_term(3, (1, 0))
        t2 = biv_term(2, (0, 1))
        ae = t1 + t2
        @test ae isa SA.AlgebraElement
        # Check it has two terms
        c = SA.coeffs(ae)
        @test length(collect(SA.keys(c))) == 2

        # Same index: should combine
        t3 = biv_term(5, (1, 0))
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

        # Rational
        tr = biv_term(3 // 4, (2, 0))
        @test coefficient(tr) == 3 // 4
        @test coefficient(zero(tr)) == 0 // 1
    end
end
