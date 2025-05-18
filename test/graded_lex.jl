# Example with Graded Lex Ordering
using Test
import StarAlgebras as SA

@testset "Graded Lex" begin
    alg = bivariate_algebra()
    o = one(alg)
    @test isone(o)
    @test SA.coeffs(o).isless == grlex
    a = SA.AlgebraElement(
        SA.SparseCoefficients(
            collect(Iterators.take(SA.object(SA.basis(alg)), 3)),
            [2, -1, 3],
            grlex,
        ),
        alg,
    )
    c = a * a
    @test c.coeffs.values == [4, -4, 12, 1, -6, 9]
    @test c.coeffs.basis_elements == [(0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0)]
end
