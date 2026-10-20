# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module TestScalarArithmetic

using Test
using SparseArrays
import StarAlgebras as SA
import MutableArithmetics as MA
import JuMP
import MathOptInterface as MOI

include("../examples/bivariate.jl")
SA.star(m::Monomial) = m

struct Scalar
    value::Rational{Int}
end

SA.@scalar_arithmetic Scalar

Base.zero(::Type{Scalar}) = Scalar(0)
Base.zero(::Scalar) = zero(Scalar)
Base.one(::Type{Scalar}) = Scalar(1)
Base.one(::Scalar) = one(Scalar)
Base.oneunit(::Type{Scalar}) = one(Scalar)
Base.iszero(a::Scalar) = iszero(a.value)
Base.broadcastable(a::Scalar) = Ref(a)
Base.:(==)(a::Scalar, b::Scalar) = a.value == b.value
Base.convert(::Type{Scalar}, a::Number) = Scalar(a)
Base.promote_rule(::Type{Scalar}, ::Type{<:Number}) = Scalar
Base.:+(a::Scalar) = a
Base.:-(a::Scalar) = Scalar(-a.value)
for op in (:+, :-, :*, :/, ://, :div)
    @eval begin
        Base.$op(a::Scalar, b::Scalar) = Scalar($op(a.value, b.value))
        Base.$op(a::Scalar, b::Number) = $op(a, Scalar(b))
        Base.$op(a::Number, b::Scalar) = $op(Scalar(a), b)
    end
end

function check_coefficients(p, expected)
    for (exp, value) in zip(((0, 0), (1, 0), (0, 1)), expected)
        i = SA.basis(parent(p))[Monomial(exp)]
        actual =
            p isa SA.Term ?
            (i == p.index ? SA.coefficient(p) : zero(SA.coefficient(p))) :
            SA.coeffs(p)[i]
        @test MA.isequal_canonical(actual, convert(typeof(actual), value))
    end
    return
end

function test_scalar(a, values)
    fixed = SA.FixedBasis(Monomial.([(0, 0), (1, 0), (0, 1)]))
    for alg in (SA.StarAlgebra(Monomial((0, 0)), fixed), bivariate_algebra())
        keys =
            [SA.basis(alg)[Monomial(exp)] for exp in ((0, 0), (1, 0), (0, 1))]
        t = SA.Term(alg, keys[2], values[2])
        containers = if SA.basis(alg) isa SA.ExplicitBasis
            (copy(values), sparsevec(1:3, copy(values), 3))
        else
            order = sortperm(keys; lt = grlex)
            (
                SA.SparseCoefficients(keys[order], values[order], grlex),
                SA.SparseCoefficients(
                    Tuple(keys[order]),
                    Tuple(values[order]),
                    grlex,
                ),
            )
        end
        for c in containers, p in (SA.AlgebraElement(c, alg), t)
            pv =
                p isa SA.Term ? [zero(values[2]), values[2], zero(values[2])] :
                values
            original = MA.copy_if_mutable(p)
            for op in (+, -, *), left in (false, true)
                args = left ? (a, p) : (p, a)
                result = @inferred op(args...)
                expected = if op === (*)
                    left ? a .* pv : pv .* a
                else
                    av = [a, zero(pv[2]), zero(pv[3])]
                    left ? op.(av, pv) : op.(pv, av)
                end
                check_coefficients(result, expected)
                @test parent(result) === alg
                @test (@inferred MA.promote_operation(op, typeof.(args)...)) ===
                      typeof(result)
                if result isa SA.AlgebraElement
                    output = zero(result)
                    @test MA.operate_to!(output, op, args...) === output
                    check_coefficients(output, expected)
                    if p isa SA.AlgebraElement && !left
                        check_coefficients(
                            MA.operate!!(op, MA.copy_if_mutable(p), a),
                            expected,
                        )
                    end
                end
            end
            for op in (MA.add_mul, MA.sub_mul), left in (false, true)
                args = left ? (a, p) : (p, a)
                product = left ? a .* pv : pv .* a
                R = MA.promote_operation(
                    op,
                    typeof(SA.algebra_element(t)),
                    typeof.(args)...,
                )
                output = MA.zero_for(R, (p,))
                check_coefficients(
                    MA.operate!(op, output, args...),
                    MA.add_sub_op(op).(product),
                )
                output = zero(output)
                f = MA.mutable_copy(output)
                @test MA.operate_to!(output, op, f, args...) === output
                check_coefficients(output, MA.add_sub_op(op).(product))
            end
            check_coefficients(p, pv)
            check_coefficients(original, pv)
        end
    end
    # Constructor dispatch must be available without loading PolyJuMP or MP.
    constant = bivariate_algebra()(a)
    @test MA.isequal_canonical(SA.coeffs(constant)[(0, 0)], a)
    @test length(collect(SA.nonzero_pairs(SA.coeffs(constant)))) == 1
    return
end

@testset "Registered scalar arithmetic" begin
    @test Base.get_extension(SA, :StarAlgebrasJuMPExt) !== nothing
    @test Base.get_extension(SA, :StarAlgebrasMathOptInterfaceExt) !== nothing
    model = JuMP.Model()
    JuMP.@variable(model, x)
    v = MOI.VariableIndex(1)
    for a in (2, 1.5, Scalar(2), x, 2x + 1, v, 2v + 1)
        @testset "$(typeof(a))" begin
            test_scalar(a, [1, 2, 3])
        end
    end
    for a in (2, x, 2x + 1)
        test_scalar(a, [x + 1, 2x, x + 3])
    end
    for a in (2, v, 2v + 1)
        test_scalar(a, [v + 1, 2v, v + 3])
    end
    test_scalar(Scalar(2), Scalar.([1, 2, 3]))
end

@testset "Scalar division and coefficient containers" begin
    alg = bivariate_algebra()
    c = SA.SparseCoefficients([(0, 0), (1, 0)], [3, 5], grlex)
    p = SA.AlgebraElement(c, alg)
    for a in (2, Scalar(2)), op in (*, /, //, div)
        result = @inferred op(p, a)
        expected = [op(3, a), op(5, a), zero(op(5, a))]
        check_coefficients(result, expected)
        output = zero(result)
        @test MA.operate_to!(output, op, p, a) === output
        check_coefficients(output, expected)
        output = convert(typeof(result), MA.mutable_copy(p))
        @test MA.operate!(op, output, a) === output
        check_coefficients(output, expected)
        coefficients = @inferred op(c, a)
        check_coefficients(SA.AlgebraElement(coefficients, alg), expected)
        if op === (*)
            coefficients = @inferred a * c
            check_coefficients(SA.AlgebraElement(coefficients, alg), expected)
        end
        check_coefficients(p, [3, 5, 0])
    end
    @test_throws MethodError p * "unregistered scalar"
    @test_throws MethodError SA.Term(alg, (1, 0), 2) + "unregistered scalar"
end

@testset "Rewritten scalar arithmetic" begin
    model = JuMP.Model()
    JuMP.@variable(model, x)
    v = MOI.VariableIndex(1)
    alg = bivariate_algebra()
    p = SA.AlgebraElement(
        SA.SparseCoefficients([(0, 0), (1, 0)], [3, 5], grlex),
        alg,
    )
    for a in (2, Scalar(2), x, 2x + 1, v, 2v + 1)
        for result in (
            MA.@rewrite(p + a * p),
            MA.@rewrite(p - p * a),
            MA.@rewrite(a + p),
            MA.@rewrite(p - a),
        )
            @test parent(result) === alg
        end
        check_coefficients(
            MA.@rewrite(p + a * p),
            [3 + a * 3, 5 + a * 5, zero(5 + a * 5)],
        )
        check_coefficients(
            MA.@rewrite(p - p * a),
            [3 - 3 * a, 5 - 5 * a, zero(5 - 5 * a)],
        )
        check_coefficients(MA.@rewrite(a + p), [a + 3, 5, 0])
        check_coefficients(MA.@rewrite(p - a), [3 - a, 5, 0])
        check_coefficients(p, [3, 5, 0])
    end
end

allocated_product(a, b) = @allocated a * b

@testset "Scalar term inference and allocations" begin
    t = SA.Term(bivariate_algebra(), (1, 0), 2)
    for a in (1.5, Scalar(2)), (left, right) in ((a, t), (t, a))
        @test SA.coefficient(@inferred left * right) == a * 2
        if VERSION >= v"1.13"
            allocated_product(left, right)
            @test allocated_product(left, right) == 0
        end
    end
end

end
