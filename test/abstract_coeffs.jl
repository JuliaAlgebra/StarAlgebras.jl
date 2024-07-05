@testset "Abstract coefficients" begin
    G = PermGroup(perm"(1,2,3)", perm"(1,2)")
    RG = StarAlgebra(G, SA.DiracBasis(G))
    fRG = let RG = RG, n = length(basis(RG))
        fb = SA.FixedBasis(basis(RG); n = n)
        StarAlgebra(SA.object(RG), fb)
    end

    h = Permutation(perm"(2,3)", G)

    α, β = let a = fRG(1), b = fRG(h)
        x = AlgebraElement(ACoeffs(coeffs(a, basis(fRG))), fRG)
        y = AlgebraElement(ACoeffs(coeffs(b, basis(fRG))), fRG)
        x, y
    end
    @test coeffs(α - β, SA.FixedBasis(basis(RG); n = 2)) == [1, -1]
    @test_throws KeyError coeffs(α - β, SA.FixedBasis(basis(RG); n = 1))
    @test SA.adjoint_coeffs(α - β, SA.FixedBasis(basis(RG); n = 2)) == [1, -1]
    @test SA.adjoint_coeffs(α - β, SA.FixedBasis(basis(RG); n = 1)) == [1]
    @test coeffs(2α) isa ACoeffs{Int}
    @test coeffs(α - β) isa ACoeffs
    @test coeffs(α - β // 3) isa ACoeffs{<:Rational}
    @test 3(α // 3) == α

    l = length(basis(RG))
    P = sum(RG(g) for g in basis(RG)) // l
    fP = AlgebraElement(ACoeffs(coeffs(P, basis(fRG))), fRG)
    @test AlgebraElement(coeffs(fP, basis(RG)), RG) == P
    @test fP * fP == fP

    fP2 = (α + β) // 2
    @test coeffs(fP2) isa ACoeffs{<:Rational}

    @test SA.key_type(coeffs(fP2)) == UInt32
    @test zero(coeffs(fP2)) == ACoeffs(coeffs(zero(P), basis(fRG)))

    P3 = 2 * sum(RG(g) for g in basis(RG) if sign(g) > 0) // l
    fP3 = AlgebraElement(ACoeffs(coeffs(P3, basis(fRG))), fRG)
    @test AlgebraElement(coeffs(fP3, basis(RG)), RG) == P3
    @test fP3 * fP3 == fP3

    PAlt = sum(sign(g) * RG(g) for g in basis(RG)) // l
    fPAlt = AlgebraElement(ACoeffs(coeffs(PAlt, basis(fRG))), fRG)
    @test AlgebraElement(coeffs(fPAlt, basis(RG)), RG) == PAlt
    @test fPAlt * fPAlt == fPAlt

    @test fP3 * fPAlt == fPAlt * fP3

    @test fP2 * fP2 == fP2

    @test fP2 * fP3 == fP3 * fP2 == fP

    P2m = (RG(1) - RG(h)) // 2
    fP2m = AlgebraElement(ACoeffs(coeffs(P2m, basis(fRG))), fRG)
    @test fP2m * fP2m == fP2m

    @test fP2m * fP3 == fP3 * fP2m == fPAlt
    @test iszero(fP2m * fP2)

    @test norm(fP2m) == norm(P2m) == norm(fP2m)
    v = coeffs(P2m, basis(fRG)) # an honest vector
    @test dot(fP2m, fP2m) == dot(coeffs(fP2m), v) == dot(v, coeffs(fP2m))

    s1, s2 = PermutationGroups.gens(G)
    @assert s1 * s2 ≠ s2 * s1
    Z = RG(1) + RG(s1)
    @test s2 * Z == RG(s2) + RG(s2 * s1)
    @test Z * s2 == RG(s2) + RG(s1 * s2)
end
