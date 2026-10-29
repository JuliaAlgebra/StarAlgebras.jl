# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

using Test, LinearAlgebra, StarAlgebras
using PermutationGroups
import StarAlgebras as SA
import MutableArithmetics as MA

if !@isdefined aug
    include("group_utils.jl")
end

if !@isdefined ACoeffs
    include("../examples/acoeffs.jl")
end

struct Wrap{T}
    x::T
end
Base.:*(a::Wrap, b::Wrap) = Wrap(a.x * b.x)

wrap(x) = Wrap(x)
unwrap(w::Wrap) = w.x

function coefficient_dot_allocations(a, b)
    dot(a, b)
    return @allocated dot(a, b)
end

@testset "Abstract coefficients" begin
    @testset "Equality with filtered iterators" begin
        a = ACoeffs([1, 0, 2])
        @test a == ACoeffs([1, 0, 2, 0])
        @test a != ACoeffs([1, 2, 0])
        @test a != ACoeffs([1, 0, 3])
        for b in (ACoeffs([1]), ACoeffs(Int[]))
            @test a != b
            @test b != a
        end
        @test ACoeffs(Int[]) == ACoeffs([0])
    end

    @testset "Coefficient dot products" begin
        sparse_coefficients(v) =
            SA.SparseCoefficients(findall(!iszero, v), filter(!iszero, v))
        for T in (Float32, Float64, ComplexF32, ComplexF64)
            a, b = T[1, 0, 2], T[3, 5, -2]
            if T <: Complex
                a += T[2im, 0, -im]
                b += T[4im, -im, 3im]
            end
            for av in (a, zero(a)),
                bv in (b, zero(b)),
                make_a in (ACoeffs, sparse_coefficients),
                make_b in (ACoeffs, sparse_coefficients)

                ac, bc = make_a(av), make_b(bv)
                @test (@inferred dot(ac, bc)) === dot(av, bv)
                @test (@inferred dot(bc, ac)) === dot(bv, av)
                @test (@inferred dot(av, bc)) === dot(av, bv)
                @test (@inferred dot(ac, bv)) === dot(av, bv)
            end
            empty_coeffs = SA.SparseCoefficients(Int[], T[])
            @test dot(empty_coeffs, empty_coeffs) === zero(T)
            @test dot(T[], empty_coeffs) === zero(T)
            @test dot(empty_coeffs, T[]) === zero(T)
        end
        # Dot recursively pairs coefficient values, including array values.
        a = SA.SparseCoefficients([1], [[1 + 2im, 3]])
        b = SA.SparseCoefficients([1], [[4 + 5im, 6]])
        @test dot(a, b) == dot(only(values(a)), only(values(b))) == 32 - 3im
    end

    @testset "Mutable coefficient dot" begin
        allocations = map((4, 256)) do n
            a, b = fill(big(2), n), fill(big(3), n)
            ac = SA.SparseCoefficients(collect(1:n), a)
            bc = SA.SparseCoefficients(collect(1:n), b)
            return map(((ac, bc), (a, bc), (ac, b))) do (left, right)
                result = @inferred dot(left, right)
                @test result == 6n
                # Only the result may be mutated, even with shared inputs.
                MA.operate!(+, result, big(1))
                @test all(==(2), a)
                @test all(==(3), b)
                @test (@inferred dot(ac, ac)) == 4n
                @test all(==(2), a)
                return coefficient_dot_allocations(left, right)
            end
        end
        # The accumulator and multiplication buffer are allocated once.
        # More matching coefficients must not add per-entry allocations.
        @test all(allocations[2] .<= allocations[1] .+ 64)
    end

    G = PermGroup(perm"(1,2,3)", perm"(1,2)")
    RG = StarAlgebra(G, SA.DiracBasis(G))
    wRG = StarAlgebra(G, SA.MappedBasis(G, wrap, unwrap))

    for basis in [RG, wRG]
        A = SA.AlgebraElement(
            SA.SparseCoefficients(collect(G), Float64.(1:6)),
            basis,
        )
        B = SA.AlgebraElement(
            SA.SparseCoefficients(collect(G)[2:2], [-1]),
            basis,
        )
        @test (A - B)^2 == A^2 - A * B - B * A + B^2
    end

    fRG = let RG = RG, n = UInt32(length(basis(RG)))
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
    @test dot(fP2m, fP2m) == adjoint(fP2m) * fP2m
    @test dot(coeffs(fP2m), coeffs(fP2m)) ==
          dot(coeffs(fP2m), v) ==
          dot(v, coeffs(fP2m))

    s1, s2 = PermutationGroups.gens(G)
    @assert s1 * s2 ≠ s2 * s1
    Z = RG(1) + RG(s1)
    @test s2 * Z == RG(s2) * Z == RG(s2) + RG(s2 * s1)
    @test Z * s2 == Z * RG(s2) == RG(s2) + RG(s1 * s2)
end

@testset "Tuple coefficients" begin
    @test SA.key_type((1, 2)) == Int
    @test SA.value_type((1, 2)) == Int
    @test SA.value_type((1.0,)) == Float64
    @test SA.value_type(NTuple{3,Rational{Int}}) == Rational{Int}
    @test SA.value_type(NTuple{2,Integer}) == Integer
end
