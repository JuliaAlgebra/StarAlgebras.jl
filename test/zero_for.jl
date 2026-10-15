# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module TestZeroFor

using Test
using SparseArrays
import StarAlgebras as SA
import MutableArithmetics as MA

@testset "Algebra zeros and reductions" begin
    alg = SA.StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0]))
    for T in (Int, BigInt),
        c in (
            T[1, 0, 2],
            sparsevec([1, 3], T[1, 2], 3),
            SA.SparseCoefficients([1, 3], T[1, 2]),
            SA.SparseCoefficients((1, 3), (T(1), T(2))),
        )

        a = SA.AlgebraElement(c, alg)
        original =
            [k => MA.copy_if_mutable(v) for (k, v) in SA.nonzero_pairs(c)]
        P = MA.promote_operation(+, typeof(a), typeof(a))
        z = @inferred MA.zero_for(P, [a, a])
        other = @inferred MA.zero_for(P, [a, a])
        @test typeof(z) === P
        @test parent(z) === alg
        @test iszero(z) && iszero(other)
        @test SA.coeffs(z) !== SA.coeffs(other)
        SA.coeffs(z)[1] = T(3)
        @test iszero(other)
        if T === BigInt && SA.coeffs(other) isa Vector
            MA.operate!(+, SA.coeffs(other)[1], 1)
            @test iszero(SA.coeffs(other)[2])
        end
        @test (@inferred sum([a])) == a
        @test (@inferred sum([a, a])) == 2a
        @test (@inferred MA.operate(sum, [a, a])) == 2a
        @test (@inferred MA.fused_map_reduce(MA.add_mul, [a, a])) == 2a
        @test sum(reshape([a, a, a, a], 2, 2); dims = 1) == fill(2a, 1, 2)
        init = zero(a)
        @test MA.operate(sum, typeof(a)[]; init) === init
        @test sum(MA.AbstractMutable[]; init) == init
        @test sum(reshape(typeof(a)[], 0, 2); dims = 1, init) ==
              fill(init, 1, 2)
        @test [k => v for (k, v) in SA.nonzero_pairs(SA.coeffs(a))] == original
        @test_throws MethodError sum(typeof(a)[])
        @test_throws MethodError MA.fused_map_reduce(MA.add_mul, typeof(a)[])
    end

    a = SA.AlgebraElement([1, 0, 2], alg)
    b = SA.AlgebraElement(SA.SparseCoefficients([1], [3]), alg)
    P = MA.promote_operation(+, typeof(b), typeof(a))
    z = @inferred MA.zero_for(P, [b], [a], [2, 3])
    @test typeof(z) === typeof(a)
    @test parent(z) === alg
    @test iszero(z)

    t = SA.Term(alg, 1, 2)
    z = @inferred MA.zero_for(typeof(zero(Int, alg)), (t,))
    @test parent(z) === alg
    @test iszero(z)
end

function test_zero_allocations(prototype, inputs::Vararg{Any,N}) where {N}
    P = MA.promote_operation(zero, typeof(prototype))
    z = @inferred MA.zero_for(P, inputs...)
    @test typeof(z) === P
    @test parent(z) === parent(prototype)
    @test iszero(z)
    @test SA.coeffs(z) !== SA.coeffs(prototype)
    # Keep both results alive so the allocation checks include their storage.
    direct, output = Ref{P}(), Ref{P}()
    direct[] = zero(prototype)
    output[] = MA.zero_for(P, inputs...)
    budget = @allocated direct[] = zero(prototype)
    @test (@allocated output[] = MA.zero_for(P, inputs...)) <= budget
    @test iszero(direct[]) && iszero(output[])
    @test SA.coeffs(direct[]) !== SA.coeffs(output[])
    return
end

@testset "Heterogeneous input inference and allocations" begin
    alg = SA.StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0]))
    for T in (Int, BigInt)
        a = SA.AlgebraElement(T[1, 0, 2], alg)
        b = SA.AlgebraElement(SA.SparseCoefficients([1, 3], T[1, 2]), alg)
        t = SA.Term(alg, 1, T(3))
        for (prototype, other) in ((a, b), (b, a))
            test_zero_allocations(prototype, [prototype])
            # Each iterable has its own concrete type. No collection of their
            # elements or extra zero should be needed while preparing the zero.
            test_zero_allocations(
                prototype,
                [other],
                (t,),
                view([prototype], :),
                Iterators.reverse([other]),
                (2, 3),
                reshape([prototype], 1, 1),
            )
            test_zero_allocations(
                prototype,
                (),
                (2, 3),
                typeof(prototype)[],
                (prototype,),
            )
        end
    end
end

@testset "Basis multiplication with a prepared destination" begin
    basis = SA.DiracBasis(["", "a", "b", "ab", "ba"])
    alg = SA.StarAlgebra("", SA.DiracMStructure(basis, *))
    f = SA.AlgebraElement(SA.SparseCoefficients(["", "a"], [1, 2]), alg)
    for (op, expected) in ((SA.__lmul, ["b", "ba"]), (SA.__rmul, ["ab", "b"]))
        output = zero(f)
        @test MA.operate_to!(output, op, "b", f) === output
        @test keys(SA.coeffs(output)) == expected
    end
end

end
