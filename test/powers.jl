# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

using Test, StarAlgebras
import StarAlgebras as SA
import MutableArithmetics as MA

if !@isdefined(biv_alg)
    include("term.jl")
end

@testset "Trivial powers preserve product promotion" begin
    b = SA.MappedBasis(0:16, map_to_chebyshev, map_from_chebyshev)
    cheby_alg = StarAlgebra(ChebyPoly(0), ChebyMStruct(b))
    for alg in (biv_alg, cheby_alg), T in (Bool, Int, BigInt)
        index = alg === biv_alg ? (1, 2) : 2
        t = Term(alg, index, one(T))
        for a in (t, SA.algebra_element(t))
            P = typeof(a * a)
            for n in (0, 1)
                result = @inferred a^n
                @test typeof(result) === P
                @test SA.algebra_element(result) ==
                      SA.algebra_element(n == 0 ? one(a) : a)
                @test parent(result) === parent(a)
                if n == 1
                    if a isa AlgebraElement
                        @test SA.coeffs(result) !== SA.coeffs(a)
                    end
                    if T === BigInt && result isa Term
                        @test coefficient(result) !== coefficient(a)
                    end
                end
            end
        end
    end
end

@testset "Trivial powers do not multiply" begin
    products = Ref(0)
    multiply = (a, b) -> begin
        products[] += 1
        return a * b
    end
    ms = SA.DiracMStructure(biv_basis, multiply)
    alg = StarAlgebra(Monomial((0, 0)), ms)
    t = Term(alg, (1, 2), 3)
    for a in (t, SA.algebra_element(t)), n in (0, 1, 2)
        products[] = 0
        result = a^n
        if n < 2
            @test products[] == 0
        else
            @test products[] > 0
        end
        @test SA.algebra_element(result) ==
              SA.algebra_element(n == 0 ? one(a) : n == 1 ? a : a * a)
    end
end
