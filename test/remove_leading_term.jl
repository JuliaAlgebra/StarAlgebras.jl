# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

@testset "Remove leading term" begin
    alg = bivariate_algebra()
    indices = [(0, 0), (1, 0), (2, 0)]
    values = BigInt[1, 2, 3]
    no_compare =
        (_, _) -> error("Removing the last coefficient must not compare keys")
    a = AlgebraElement(
        SA.SparseCoefficients(copy(indices), copy(values), no_compare),
        alg,
    )
    saved = Term(alg, indices[end], values[end])

    b = remove_leading_term(a)
    @test parent(b) === alg
    @test SA.coeffs(b).isless === no_compare
    @test keys(SA.coeffs(b)) == indices[1:2]
    @test keys(SA.coeffs(a)) == indices
    @test MA.promote_operation(remove_leading_term, typeof(a)) === typeof(b)
    @test keys(SA.coeffs(MA.operate(remove_leading_term, a))) == indices[1:2]
    for remaining in 2:-1:0
        @test MA.operate!(remove_leading_term, a) === a
        @test keys(SA.coeffs(a)) == indices[1:remaining]
        @test SA.coefficient(saved) === values[end]
        @test SA.coefficient(saved) == 3
    end
    @test MA.operate!(remove_leading_term, a) === a
    @test iszero(a)
    @test iszero(MA.operate(remove_leading_term, saved))
    @test MA.promote_operation(remove_leading_term, typeof(saved)) ===
          typeof(remove_leading_term(saved))

    finite_alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0, 4.0]))
    for c in (BigInt[1, 0, 3, 0], sparsevec([1, 3], BigInt[1, 3], 4))
        a = AlgebraElement(c, finite_alg)
        saved_coefficient = c[3]
        b = remove_leading_term(a)
        @test parent(b) === finite_alg
        @test SA.coeffs(b) == [1, 0, 0, 0]
        @test SA.coeffs(a) == [1, 0, 3, 0]
        @test MA.operate!(remove_leading_term, a) === a
        @test SA.coeffs(a) == [1, 0, 0, 0]
        @test saved_coefficient == 3
        MA.operate!(remove_leading_term, a)
        @test iszero(a)
        @test MA.operate!(remove_leading_term, a) === a
        @test length(SA.coeffs(a)) == 4
        if c isa SparseVector
            @test isempty(SparseArrays.nonzeroinds(c))
            @test isempty(SparseArrays.nonzeros(c))
        end
    end
end
