struct ACoeffs{T} <: SA.AbstractCoefficients{UInt32,T}
    vals::Vector{T}
end

# convienience:
ACoeffs(v::AbstractVector) = ACoeffs{valtype(v)}(v)

# AbstractCoefficients API

## Basic API
Base.keys(ac::ACoeffs) = (k for (k, v) in pairs(ac.vals) if !iszero(v))
Base.values(ac::ACoeffs) = (v for v in ac.vals if !iszero(v))
SA.canonical(ac::ACoeffs) = ac
function SA.star(b::SA.AbstractBasis, ac::ACoeffs)
    return ACoeffs([ac.vals[star(b, k)] for k in keys(ac)])
end

# for preallocation in arithmetic
function Base.similar(ac::ACoeffs, ::Type{T}) where {T}
    # it is needed to zero as uninitialized values may not be cleared by
    # MA.operate!(zero, ac)
    return ACoeffs(zeros(T, size(ac.vals)))
end

# the default arithmetic implementation uses this access
Base.getindex(ac::ACoeffs, idx) = ac.vals[idx]
Base.setindex!(ac::ACoeffs, val, idx) = ac.vals[idx] = val

# a very basic * using ms.structure;
# TODO: this is so generic that it should be included in SA
# but then it hijacks the more refined version
function MA.operate!(
    ms::SA.UnsafeAddMul,
    res::ACoeffs,
    v::SA.AbstractCoefficients,
    w::SA.AbstractCoefficients,
)
    for (kv, a) in SA.nonzero_pairs(v)
        for (kw, b) in SA.nonzero_pairs(w)
            c = ms.structure(kv, kw)
            for (k, v) in SA.nonzero_pairs(c)
                res[ms.structure[k]] += v * a * b
            end
        end
    end
    return res
end

@testset "Abstract coefficients" begin
    G = PermGroup(perm"(1,2,3)", perm"(1,2)")
    RG = StarAlgebra(G, SA.DiracBasis{UInt8}(G))
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
end
