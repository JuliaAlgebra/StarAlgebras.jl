# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

# Mutable wrapper for testing MA operations
mutable struct MutableInt <: Number
    value::Int
end
Base.zero(::MutableInt) = MutableInt(0)
Base.zero(::Type{MutableInt}) = MutableInt(0)
Base.one(::MutableInt) = MutableInt(1)
Base.one(::Type{MutableInt}) = MutableInt(1)
Base.iszero(x::MutableInt) = iszero(x.value)
Base.:(==)(a::MutableInt, b::MutableInt) = a.value == b.value
Base.:*(a::MutableInt, b::MutableInt) = MutableInt(a.value * b.value)
MA.mutability(::Type{MutableInt}) = MA.IsMutable()
function MA.operate_to!(
    res::MutableInt,
    ::typeof(*),
    a::MutableInt,
    b::MutableInt,
)
    res.value = a.value * b.value
    return res
end
function MA.operate!(::typeof(*), a::MutableInt, b::MutableInt)
    a.value *= b.value
    return a
end
function MA.operate!(::typeof(one), a::MutableInt)
    a.value = 1
    return a
end
function MA.copy_if_mutable(a::MutableInt)
    return MutableInt(a.value)
end
function MA.mutable_copy(a::MutableInt)
    return MutableInt(a.value)
end

@testset "Term" begin
    @testset "Accessors" begin
        t = Term(3.0, Monomial((1, 2)))
        @test coefficient(t) == 3.0
        @test basis_element(t) == Monomial((1, 2))
    end

    @testset "iszero / zero" begin
        t = Term(3, Monomial((1, 0)))
        @test !iszero(t)
        z = zero(t)
        @test iszero(z)
        @test coefficient(z) == 0
        @test basis_element(z) == Monomial((1, 0))

        @test iszero(Term(0.0, Monomial((0, 0))))
        @test !iszero(Term(1, Monomial((0, 0))))
    end

    @testset "copy / mutable_copy" begin
        m = Monomial((2, 3))
        t = Term(5, m)
        t2 = copy(t)
        @test coefficient(t2) == 5
        @test basis_element(t2) == m

        t3 = MA.mutable_copy(t)
        @test coefficient(t3) == coefficient(t)
        @test basis_element(t3) == basis_element(t)

        # Test that copy of mutable fields produces independent copies
        a = MutableInt(7)
        b = MutableInt(3)
        t_mut = Term(a, b)
        t_copy = copy(t_mut)
        @test coefficient(t_copy) == MutableInt(7)
        @test basis_element(t_copy) == MutableInt(3)
        # Mutating the copy should not affect the original
        t_copy.coefficient.value = 99
        @test coefficient(t_mut) == MutableInt(7)
    end

    @testset "mutability" begin
        # Int and Monomial (immutable struct) are both not mutable
        @test MA.mutability(Term{Int,Monomial}) isa MA.IsNotMutable
        # MutableInt on both sides gives mutable
        @test MA.mutability(Term{MutableInt,MutableInt}) isa MA.IsMutable
        # Mixed: one immutable makes the term not mutable
        @test MA.mutability(Term{MutableInt,Int}) isa MA.IsNotMutable
        @test MA.mutability(Term{Int,MutableInt}) isa MA.IsNotMutable
    end

    # MA.operate_to!(*, ...), MA.operate!(*, ...) and MA.operate!(one, ...)
    # are defined by downstream packages (e.g., MultivariatePolynomials)
    # and tested there.

    @testset "one / isone" begin
        t = Term(3, Monomial((1, 0)))
        @test !isone(t)
        t1 = Term(1, Monomial((0, 0)))
        @test isone(t1)
        @test one(t1) == t1
    end

    @testset "== / isequal" begin
        t1 = Term(3.0, Monomial((1, 0)))
        t2 = Term(3.0, Monomial((1, 0)))
        t3 = Term(4.0, Monomial((1, 0)))
        t4 = Term(3.0, Monomial((0, 1)))
        @test t1 == t2
        @test t1 != t3
        @test t1 != t4
        @test isequal(t1, t2)
        @test !isequal(t1, t3)
        # zero terms are equal regardless of basis element
        @test Term(0.0, Monomial((1, 0))) == Term(0.0, Monomial((0, 1)))
    end

    @testset "hash" begin
        t1 = Term(3, Monomial((1, 0)))
        t2 = Term(3, Monomial((1, 0)))
        @test hash(t1) == hash(t2)
    end

    @testset "convert / promote_rule" begin
        t = Term(3, Monomial((1, 0)))
        t2 = convert(Term{Float64,Monomial}, t)
        @test coefficient(t2) == 3.0
        @test basis_element(t2) == Monomial((1, 0))
        @test promote_type(Term{Int,Monomial}, Term{Float64,Monomial}) ==
              Term{Float64,Monomial}
    end

    @testset "broadcastable / ndims" begin
        t = Term(2, Monomial((1, 0)))
        @test ndims(t) == 0
        @test ndims(typeof(t)) == 0
        @test Base.broadcastable(t) isa Ref
    end

    @testset "various coefficient types" begin
        # Float64 coefficient
        t = Term(2.5, Monomial((1, 1)))
        @test coefficient(t) == 2.5
        @test !iszero(t)

        # Complex coefficient
        tc = Term(1 + 2im, Monomial((0, 1)))
        @test coefficient(tc) == 1 + 2im
        @test !iszero(tc)
        @test iszero(zero(tc))
        @test coefficient(zero(tc)) == 0 + 0im

        # Rational coefficient
        tr = Term(3 // 4, Monomial((2, 0)))
        @test coefficient(tr) == 3 // 4
        @test coefficient(zero(tr)) == 0 // 1
    end
end
