module TestBasisConversion

using Test
import StarAlgebras as SA
import MutableArithmetics as MA

struct ScaledBasis{W} <: SA.ExplicitBasis{W,Int}
    weight::W
end

Base.length(::ScaledBasis) = 2
Base.getindex(b::ScaledBasis, i::Int) = b.weight * i

function SA.promote_coefficient_type(
    op::Union{typeof(SA.coeffs),typeof(SA.adjoint_coeffs)},
    ::Type{T},
    ::Type{ScaledBasis{S}},
    ::Type{ScaledBasis{D}},
) where {T,S,D}
    W = MA.promote_operation(/, S, D)
    if op === SA.adjoint_coeffs
        W = MA.promote_operation(conj, W)
    end
    return MA.promote_operation(*, W, T)
end

function SA.coeffs!(res, cfs, source::ScaledBasis, target::ScaledBasis)
    res .= (source.weight / target.weight) .* cfs
    return res
end

function SA.adjoint_coeffs!(res, cfs, source::ScaledBasis, target::ScaledBasis)
    res .= conj(source.weight / target.weight) .* cfs
    return res
end

function test_promotion(
    op,
    ::Type{T},
    source::S,
    target::D,
    expected,
) where {T,S,D}
    @test (@inferred SA.promote_coefficient_type(op, T, S, D)) === expected
    SA.promote_coefficient_type(op, T, S, D)
    # MA's mixed BigInt arithmetic promotion constructs sample scalar values.
    if isbitstype(T)
        @test (@allocated SA.promote_coefficient_type(op, T, S, D)) == 0
    end
    return
end

@testset "Reindexing preserves coefficient types" begin
    source = SA.FixedBasis([1, 2])
    target = SA.FixedBasis([2, 1, 3])
    for T in (Bool, Int, BigInt)
        cfs = T[1, 0]
        forward = SA.coeffs(cfs, source, target)
        @test forward == T[0, 1, 0]
        @test SA.value_type(forward) === T
        backward = SA.adjoint_coeffs(forward, source, target)
        @test backward == cfs
        @test SA.value_type(backward) === T
        for op in (SA.coeffs, SA.adjoint_coeffs)
            test_promotion(op, T, source, target, T)
            @test op(cfs, source, source) === cfs
            @test op(cfs, source, SA.FixedBasis([1, 2])) == cfs
        end
    end
end

@testset "Weighted conversion and its adjoint" begin
    source = ScaledBasis(1 // 1)
    for weight in (2 // 1, complex(1 // 1, 1 // 1)), T in (Int, BigInt, Float64)
        target = ScaledBasis(weight)
        cfs = T[1, 3]
        factor = source.weight / target.weight
        expected = factor .* cfs
        forward = SA.coeffs(cfs, source, target)
        @test forward == expected
        @test SA.value_type(forward) === eltype(expected)
        test_promotion(SA.coeffs, T, source, target, eltype(expected))
        expected_adjoint = conj(factor) .* cfs
        backward = SA.adjoint_coeffs(cfs, source, target)
        @test backward == expected_adjoint
        @test SA.value_type(backward) === eltype(expected_adjoint)
        test_promotion(
            SA.adjoint_coeffs,
            T,
            source,
            target,
            eltype(expected_adjoint),
        )
        @test cfs == T[1, 3]
        for op in (SA.coeffs, SA.adjoint_coeffs)
            empty_result = op(zeros(T, 2), source, target)
            @test iszero(empty_result)
            @test SA.value_type(empty_result) === eltype(expected)
        end
    end
end

end # module
