module TestIsapprox

using Test
import LinearAlgebra
import StarAlgebras as SA

@testset "Array isapprox" begin
    alg = SA.StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0]))
    for T in (Int, Rational{Int}, Float32, Float64, ComplexF64),
        make in (c -> SA.Term(alg, 2, c), c -> SA.AlgebraElement(T[0, c], alg))

        p = make(T(2))
        q = make(T(3))
        @test [p] ≈ [p]
        @test [p p; p p] ≈ [p p; p p]
        @test !([p] ≈ [q])
        @test isapprox([p], [q]; atol = 2)
        @test !isapprox([p], [q]; atol = 0.5)
        @test isapprox([p], [q]; rtol = 0.5)
        @test !isapprox([p], [q]; rtol = 0.1)
    end
    for T in (Float32, Float64),
        make in (c -> SA.Term(alg, 2, c), c -> SA.AlgebraElement(T[0, c], alg))

        p = make(one(T))
        q = make(one(T) + eps(T))
        @test [p] ≈ [q]
        @test !isapprox([p], [q]; rtol = 0)
    end
    t = SA.Term(alg, 2, 3 + 4im)
    @test LinearAlgebra.norm(t) == 5
    @test LinearAlgebra.norm(t, 1) == 5
    @test LinearAlgebra.norm(t, Inf) == 5
end

end
