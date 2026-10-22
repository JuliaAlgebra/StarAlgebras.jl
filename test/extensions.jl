# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module TestExtensions

using Test
using SparseArrays
import JuMP
import MathOptInterface as MOI
import StarAlgebras as SA

include("../examples/bivariate.jl")
SA.star(m::Monomial) = m

function test_value()
    mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    model = JuMP.direct_model(mock)
    JuMP.@variable(model, a)
    JuMP.@variable(model, b)
    MOI.Utilities.mock_optimize!(mock, [2.0, 3.0])
    value = Base.Fix1(getindex, Dict(a => 2.0, b => 3.0))
    fixed = SA.FixedBasis(Monomial.([(0, 0), (1, 0)]))
    for alg in (SA.StarAlgebra(Monomial((0, 0)), fixed), bivariate_algebra())
        keys = [SA.basis(alg)[Monomial(exp)] for exp in ((0, 0), (1, 0))]
        for (cs, expected) in (
            ([a, b], [2.0, 3.0]),
            ([2a + 1, -b + 2], [5.0, -1.0]),
            ([a * b + 1, b^2 - 1], [7.0, 8.0]),
            ([3a - 2b, b - 3], [0.0, 0.0]),
        )
            t = SA.Term(alg, keys[2], cs[2])
            for result in
                ((@inferred JuMP.value(value, t)), (@inferred JuMP.value(t)))
                @test result isa SA.Term{Float64}
                @test parent(result) === alg
                @test result.index == t.index
                @test SA.coefficient(result) == expected[2]
            end
            containers = if SA.basis(alg) isa SA.ExplicitBasis
                (cs, sparsevec(1:2, cs, 2))
            else
                (
                    SA.SparseCoefficients(keys, cs, grlex),
                    SA.SparseCoefficients(Tuple(keys), Tuple(cs), grlex),
                )
            end
            for c in containers
                p = SA.AlgebraElement(c, alg)
                for result in (
                    (@inferred JuMP.value(value, p)),
                    (@inferred JuMP.value(p)),
                )
                    @test result isa SA.AlgebraElement{Float64}
                    @test parent(result) === alg
                    @test [SA.coeffs(result)[key] for key in keys] == expected
                    @test SA.coeffs(result) !== c
                end
                @test SA.coeffs(p) === c
            end
        end
    end
    return
end

function test_printing()
    alg = SA.StarAlgebra(1, SA.FixedBasis([1, 2]))
    p = SA.AlgebraElement([1, 2], alg)
    @test JuMP.function_string(MIME"text/plain"(), p) == "1·1 + 2·2"
    @test JuMP.function_string(MIME"text/latex"(), p) ==
          "1 \\cdot 1 + 2 \\cdot 2"
    return
end

function test_without_indices()
    alg = bivariate_algebra()
    t = SA.Term(alg, (1, 0), 2)
    p = SA.algebra_element(t)
    basis = SA.basis(alg)
    fixed = SA.FixedBasis([Monomial((0, 0))])
    for op in (MOI.Utilities.map_indices, MOI.Utilities.substitute_variables)
        for x in (t, p, basis, fixed),
            container in (x, (x, x), [x], [(x, x)], typeof(x)[])
            # The callback must not be called for objects without MOI indices.
            @test (@inferred op(error, container)) === container
        end
        for x in ((), Tuple{}[], (t, p, basis), [(t, p, basis)])
            @test (@inferred op(error, x)) === x
        end
        # An MOI coefficient is not index-free, unlike a numeric coefficient.
        t_index = SA.Term(alg, (1, 0), MOI.VariableIndex(1))
        @test !applicable(op, error, t_index)
        @test !applicable(op, error, SA.algebra_element(t_index))
    end
    return
end

@testset "JuMP and MOI extensions" begin
    test_value()
    test_printing()
    test_without_indices()
end

end
