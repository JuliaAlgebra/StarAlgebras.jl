# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module TestMutableArithmetics

using Test
using SparseArrays
using LinearAlgebra
using StarAlgebras
import StarAlgebras as SA
import MutableArithmetics as MA

include(joinpath(@__DIR__, "..", "examples", "acoeffs.jl"))

@testset "In-place one" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([2.0, 1.0, 3.0]))
    for T in (Int, BigInt), empty in (false, true)
        vals = empty ? zeros(T, 3) : T[2, 4, 3]
        for c in (
            copy(vals),
            sparse(vals),
            SA.SparseCoefficients(
                findall(!iszero, vals),
                filter(!iszero, vals),
            ),
        )
            a = AlgebraElement(c, alg)
            @test (@inferred MA.operate!(one, a)) === a
            @test parent(a) === alg
            @test coeffs(a) === c
            @test isone(a)
            @test [i => v for (i, v) in SA.nonzero_pairs(c)] == [2 => one(T)]
        end
    end
    alg = StarAlgebra(1.0, SA.FixedBasis([2.0, 3.0]))
    a = AlgebraElement([2, 3], alg)
    @test_throws KeyError MA.operate!(one, a)
    @test coeffs(a) == [2, 3]
end

@testset "In-place multiplication" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 4.0]))
    for T in (Int, BigInt),
        c in (
            T[1, 2, 0],
            sparsevec([1, 2], T[1, 2], 3),
            SA.SparseCoefficients([1, 2], T[1, 2]),
        )

        a = AlgebraElement(c, alg)
        b = copy(a)
        coeffs(b)[1], coeffs(b)[2] = T(3), T(-4)
        originals = MA.mutable_copy.((a, b))
        for alias in (false, true)
            output = MA.mutable_copy(a)
            storage = coeffs(output)
            old_values = collect(values(storage))
            right = alias ? output : b
            expected = alias ? T[1, 4, 4] : T[3, 2, -8]
            @test (@inferred MA.mul!!(output, right)) === output
            @test [coeffs(output)[i] for i in 1:3] == expected
            @test coeffs(output) === storage
            @test parent(output) === alg
            @test old_values == collect(values(c))
            @test (a, b) == originals
        end
        output = copy(a)
        @test MA.mul!!(output, zero(b)) === output
        @test iszero(output)
        @test MA.mul!!(output, b) === output
        @test iszero(output)
    end
end

@testset "Mutable copies and operation routing" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0, 4.0]))
    for c in (
        BigInt[1, 0, 3, 0],
        sparsevec([1, 3], BigInt[1, 3], 4),
        SA.SparseCoefficients([1, 3], BigInt[1, 3]),
    )
        a = AlgebraElement(c, alg)
        @test (@inferred MA.promote_operation(zero, typeof(a))) ===
              typeof(zero(a))
        original = deepcopy(a)
        b = MA.copy_if_mutable(a)
        @test b !== a
        @test parent(b) === parent(a)
        @test SA.coeffs(b) !== c
        @test SA.coeffs(b)[1] !== c[1]
        MA.operate!(+, SA.coeffs(b)[1], 1)
        @test SA.coeffs(b)[1] == 2
        @test c[1] == 1

        @test MA.operate!!(remove_leading_term, b) === b
        @test iszero(SA.coeffs(b)[3])
        @test SA.coeffs(b)[1] == 2
        @test MA.operate!!(zero, b) === b
        @test iszero(b)
        @test a == original

        @test (@inferred MA.add!!(b, a)) === b
        @test b == original
        @test (@inferred MA.add!!(b, b)) === b
        @test b == 2a
        @test iszero(MA.add!!(b, -b))

        @test (@inferred MA.operate!!(MA.sub_mul, a, 2, a)) === a
        @test a == -original
    end
end

@testset "Fused scalar products" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0, 4.0, 6.0, 9.0]))
    for T in (Int, BigInt), storage in (:dense, :sparse, :coefficients)
        make_element = function (v)
            c = if storage === :dense
                T.(v)
            elseif storage === :sparse
                sparse(T.(v))
            else
                keys = findall(!iszero, v)
                SA.SparseCoefficients(keys, T.(v[keys]))
            end
            return AlgebraElement(c, alg)
        end
        f = make_element([1, 0, -1, 0, 0, 0])
        a = make_element([1, 2, 0, 0, 0, 0])
        b = make_element([2, 0, 3, 0, 0, 0])
        for op in (MA.add_mul, MA.sub_mul),
            (factors, product) in (
                ((2, 3), make_element([6, 0, 0, 0, 0, 0])),
                ((2I, 3, -I), make_element([-6, 0, 0, 0, 0, 0])),
                ((a, b, 2), make_element([4, 8, 6, 0, 12, 0])),
            )

            output = MA.mutable_copy(f)
            originals = MA.copy_if_mutable.(factors)
            @test (@inferred MA.operate!(op, output, factors...)) === output
            @test output == MA.add_sub_op(op)(f, product)
            @test factors == originals
        end
        for op in (MA.add_mul, MA.sub_mul),
            (left, right) in (
                (a, b),
                (2, b),
                (b, 2),
                (2, SA.Term(alg, 2, T(3))),
                (SA.Term(alg, 2, T(3)), 2),
            )

            expected = op(
                f,
                left isa Term ?
                make_element([
                    i == left.index ? coefficient(left) : zero(T) for i in 1:6
                ]) : left,
                right isa Term ?
                make_element([
                    i == right.index ? coefficient(right) : zero(T) for i in 1:6
                ]) : right,
            )
            originals = MA.copy_if_mutable.((f, left, right))
            output = zero(f)
            @test (@inferred MA.operate_to!(output, op, f, left, right)) ===
                  output
            @test output == expected
            @test (f, left, right) == originals
            output = MA.mutable_copy(f)
            @test (@inferred MA.operate!!(op, output, left, right)) === output
            @test output == expected
        end
        for op in (MA.add_mul, MA.sub_mul), alias in (:left, :right, :both)
            output = MA.mutable_copy(f)
            left = alias === :right ? a : output
            right = alias === :left ? b : output
            expected = op(output, left, right)
            @test (@inferred MA.operate!!(op, output, left, right)) === output
            @test output == expected
        end
        output = MA.mutable_copy(b)
        for (left, right) in ((output, a), (a, output))
            @test_throws ArgumentError MA.operate_to!(
                output,
                MA.add_mul,
                f,
                left,
                right,
            )
        end
        @test output == b
        @test_throws InexactError MA.operate!(MA.add_mul, output, 0.5, a)
        # A failed coefficient write may have changed earlier output entries.
        @test b == make_element([2, 0, 3, 0, 0, 0])
        # Distinguish the scalar from this example's Float64 basis elements.
        for op in (MA.add_mul, MA.sub_mul),
            (left, right) in (
                (0.5f0, a),
                (a, 0.5f0),
                (0.5f0, SA.Term(alg, 2, T(3))),
                (SA.Term(alg, 2, T(3)), 0.5f0),
            )

            g = left isa Number ? right : left
            g = g isa Term ? make_element([0, 3, 0, 0, 0, 0]) : g
            expected = op(f, 0.5f0, g)
            @test (@inferred MA.promote_operation(
                op,
                typeof(f),
                typeof(left),
                typeof(right),
            )) === typeof(expected)
            output = zero(expected)
            @test (@inferred MA.operate_to!(output, op, f, left, right)) ===
                  output
            @test output == expected
            output = MA.mutable_copy(convert(typeof(expected), f))
            @test (@inferred MA.operate!!(op, output, left, right)) === output
            @test output == expected
        end
        # This algebra has no type-only zero; array products need its values.
        @test_throws MethodError zero(typeof(a))
        @test (@inferred MA.operate(*, [a b; b a], [a, b])) ==
              [a * a + b * b, b * a + a * b]
        @test (@inferred MA.fused_map_reduce(MA.add_mul, [a, b], [a, b])) ==
              a * a + b * b
        source = make_element([0, 0, 0, 0, 0, 1])
        output = zero(source)
        MA.operate_to!(output, MA.add_mul, source, a, b)
        coeffs(output)[6] = MA.add!!(coeffs(output)[6], 1)
        @test coeffs(source)[6] == 1
        output = MA.mutable_copy(f)
        @test MA.operate_to!(output, MA.add_mul, output, a, b) === output
        @test output == f + a * b
    end
end

# Coefficients are only scalars if their type is registered
SA.@scalar_arithmetic Matrix{Int}

@testset "Scalar coefficient order and ownership" begin
    alg = StarAlgebra("", SA.DiracBasis(["", "a", "b", "ab", "ba"]))
    A, B = [1 2; 0 1], [1 0; 3 1]
    g = AlgebraElement(SA.SparseCoefficients(["a"], [B]), alg)
    original = MA.mutable_copy(g)
    for op in (MA.add_mul, MA.sub_mul), left in (true, false)
        output = zero(g)
        a, b = left ? (A, g) : (g, A)
        @test (@inferred MA.operate!!(op, output, a, b)) === output
        expected = MA.add_sub_op(op)(left ? A * B : B * A)
        @test coeffs(output)["a"] == expected
        @test g == original
        @test A == [1 2; 0 1]
    end
    g = AlgebraElement(SA.SparseCoefficients(["a"], [big(2)]), alg)
    scalar = coeffs(g)["a"]
    @test MA.operate!(MA.add_mul, g, scalar, g) === g
    @test coeffs(g)["a"] == 6
    @test scalar == 2

    a = AlgebraElement(SA.SparseCoefficients(["a"], [A]), alg)
    b = AlgebraElement(SA.SparseCoefficients(["b"], [B]), alg)
    output = zero(a)
    @test (@inferred MA.operate!!(MA.add_mul, output, a, b)) === output
    @test coeffs(output)["ab"] == A * B
    output = copy(a)
    @test (@inferred MA.mul!!(output, b)) === output
    @test coeffs(output)["ab"] == A * B
    @test A == [1 2; 0 1]
    @test B == [1 0; 3 1]

    a = AlgebraElement(SA.SparseCoefficients(["a"], [2]), alg)
    b = AlgebraElement(SA.SparseCoefficients(["b"], [3]), alg)
    for (left, right, key) in ((a, b, "ab"), (b, a, "ba"))
        output = zero(a)
        @test (@inferred MA.operate!(MA.add_mul, output, left, right)) ===
              output
        @test keys(coeffs(output)) == [key]
        @test values(coeffs(output)) == [6]
        output = copy(left)
        @test (@inferred MA.mul!!(output, right)) === output
        @test keys(coeffs(output)) == [key]
        @test values(coeffs(output)) == [6]
    end
end

@testset "Algebra-valued scalar coefficients" begin
    inner = StarAlgebra("", SA.DiracBasis(["", "a", "b", "ab", "ba"]))
    a = AlgebraElement(SA.SparseCoefficients(["a"], [2]), inner)
    b = AlgebraElement(SA.SparseCoefficients(["b"], [3]), inner)
    outer = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 4.0]))
    p = AlgebraElement(SA.SparseCoefficients([2], [b]), outer)
    originals = MA.mutable_copy.((a, b, p))
    for (left, right, key) in ((a, p, "ab"), (p, a, "ba"))
        result = @inferred left * right
        @test (@inferred MA.promote_operation(
            *,
            typeof(left),
            typeof(right),
        )) === typeof(result)
        @test parent(result) === outer
        @test keys(coeffs(result)) == [2]
        expected = AlgebraElement(SA.SparseCoefficients([key], [6]), inner)
        @test coeffs(result)[2] == expected
        output = zero(p)
        @test (@inferred MA.operate_to!(output, *, left, right)) === output
        @test output == result
        @test (a, b, p) == originals
        MA.operate!(MA.add_mul, coeffs(result)[2], 2, 3)
        @test (a, b, p) == originals
    end
end

@testset "Unsorted scaled accumulation" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0]))
    g = AlgebraElement(SA.SparseCoefficients([3, 1, 3], [1, 2, 3], >), alg)
    f = AlgebraElement(SA.SparseCoefficients([1, 3], [1, 2]), alg)
    @test (@inferred MA.operate!(MA.sub_mul, f, 2, g)) === f
    @test keys(coeffs(f)) == [1, 3]
    @test values(coeffs(f)) == [-3, -6]
    @test keys(coeffs(g)) == [3, 1, 3]
    @test values(coeffs(g)) == [1, 2, 3]
    @test MA.operate!(MA.add_mul, g, g, 2) === g
    @test keys(coeffs(g)) == [3, 1]
    @test values(coeffs(g)) == [12, 6]
end

function test_scalar_product_allocations(::Type{T}) where {T}
    alg = StarAlgebra(1, SA.FixedBasis(collect(1:100)))
    for c in (
        ones(T, 100),
        sparse(ones(T, 100)),
        SA.SparseCoefficients(collect(1:100), ones(T, 100)),
    )
        f, g = AlgebraElement(copy(c), alg), AlgebraElement(copy(c), alg)
        for op in (MA.add_mul, MA.sub_mul)
            @test (@inferred MA.operate!(op, f, 2, g)) === f
            @test (@allocated MA.operate!(op, f, 2, g)) == 0
            @test all(isone, values(c))
        end
    end
    return
end

@testset "Scaled accumulation allocations" begin
    test_scalar_product_allocations(Int)
    test_scalar_product_allocations(Float64)
end

function test_addition_allocations(::Type{T}, n) where {T}
    alg = StarAlgebra(1, SA.FixedBasis(collect(1:max(1, n))))
    c = SA.SparseCoefficients(collect(1:n), ones(T, n))
    d = SA.SparseCoefficients(collect(1:n), -ones(T, n))
    f, g = AlgebraElement(c, alg), AlgebraElement(d, alg)
    sizehint!(keys(c), 2n)
    sizehint!(values(c), 2n)
    @test (@inferred MA.add!!(f, g)) === f
    resize!(keys(c), n)
    resize!(values(c), n)
    copyto!(keys(c), 1:n)
    fill!(values(c), one(T))
    @test (@allocated MA.add!!(f, g)) == 0
    @test iszero(f)
    @test values(d) == -ones(T, n)
    return
end

@testset "Allocation-free algebra element addition" begin
    for T in (Int, Float64), n in (0, 1, 100)
        test_addition_allocations(T, n)
    end
end

@testset "Zero promotion with tuple storage" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0, 4.0]))
    a = AlgebraElement(SA.SparseCoefficients((1, 3), (big(1), big(3))), alg)
    z = zero(a)
    @test (@inferred MA.promote_operation(zero, typeof(a))) === typeof(z)
    @test typeof(z) !== typeof(a)
    @test parent(z) === alg

    r = @inferred MA.operate!!(zero, a)
    @test iszero(r)
    @test typeof(r) === typeof(z)
    @test parent(r) === alg
    @test keys(SA.coeffs(a)) == (1, 3)
    @test values(SA.coeffs(a)) == (1, 3)
end

@testset "Coefficient storage mutation routing" begin
    v = sparsevec([1], BigInt[2], 3)
    @test MA.copy_if_mutable(v) === v

    c = SA.SparseCoefficients([1, 3], BigInt[2, 3])
    b = MA.copy_if_mutable(c)
    @test b !== c
    MA.operate!(+, b[1], 1)
    @test c[1] == 2
    result = MA.operate!!(+, b, c)
    @test result === b
    @test result[1] == 5
    @test c[1] == 2
    @test MA.operate!!(remove_leading_term, b) === b
    @test iszero(b[3])
    @test MA.operate!!(zero, b) === b
    @test iszero(b)
end

# Whether `a` stores one of the mutable values of `b`
_shares_values(a, b) = any(v -> any(Base.Fix1(===, v), values(b)), values(a))

@testset "Sparse coefficient merge" begin
    c = SA.SparseCoefficients([1, 3], [2, 4])
    pairs = SA._CoefficientPairs(c)
    @test setindex!(pairs, 2 => 7, 1) === pairs
    @test keys(c) == [2, 3]
    @test values(c) == [7, 4]

    for T in (Int, BigInt), alias in (:none, :left, :right)
        x = SA.SparseCoefficients([1, 3, 5], T[2, 3, 4])
        y = SA.SparseCoefficients([2, 3, 6], T[5, -3, 7])
        out = alias === :left ? x : alias === :right ? y : zero(x)
        @test MA.operate_to!(out, +, x, y) === out
        @test keys(out) == [1, 2, 5, 6]
        @test values(out) == [2, 5, 4, 7]
        alias === :left || @test values(x) == [2, 3, 4]
        alias === :right || @test values(y) == [5, -3, 7]
        if T == BigInt
            # Modifying `out` in place must not modify the inputs.
            alias === :left || @test !_shares_values(out, x)
            alias === :right || @test !_shares_values(out, y)
        end
        @test MA.operate!!(+, out, out) === out
        @test values(out) == [4, 10, 8, 14]
        @test MA.operate!!(+, out, -out) === out
        @test isempty(keys(out))
        @test isempty(values(out))
    end

    x = SA.SparseCoefficients([3, 1, 3], [1, 2, -1])
    y = SA.SparseCoefficients([2, 1], [3, -2], >)
    @test MA.operate!!(+, x, y) === x
    @test keys(x) == [2]
    @test values(x) == [3]
    @test keys(y) == [2, 1]
    @test values(y) == [3, -2]

    # Unsorted storage takes another path, the result must also own its values
    for alias in (:none, :left, :right)
        x = SA.SparseCoefficients([3, 1], big.([1, 2]))
        y = SA.SparseCoefficients([2], big.([3]))
        out = alias === :left ? x : alias === :right ? y : zero(x)
        @test MA.operate_to!(out, +, x, y) === out
        @test collect(SA.nonzero_pairs(out)) == [1 => 2, 2 => 3, 3 => 1]
        alias === :left || @test !_shares_values(out, x)
        alias === :right || @test !_shares_values(out, y)
    end
    out = SA.SparseCoefficients(Int[], BigInt[])
    y = SA.SparseCoefficients([2], big.([3]))
    @test MA.operate!(SA.UnsafeAdd(), out, y) === out
    @test !_shares_values(out, y)

    x = SA.SparseCoefficients([1], [1])
    y = SA.SparseCoefficients([1, 2], [0.5, 2.5])
    out = MA.operate!!(+, x, y)
    @test out !== x
    @test values(out) == [1.5, 2.5]
    @test values(x) == [1]
end

@testset "In-place negation" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([2.0, 1.0, 3.0]))
    for T in (Int, BigInt), empty in (false, true)
        vals = empty ? zeros(T, 3) : T[2, 0, 3]
        for c in (
            copy(vals),
            sparse(vals),
            SA.SparseCoefficients(
                findall(!iszero, vals),
                filter(!iszero, vals),
            ),
        )
            a = AlgebraElement(c, alg)
            @test MA.promote_operation(-, typeof(a)) == typeof(a)
            pairs = [i => v for (i, v) in SA.nonzero_pairs(c)]
            expected = [i => -v for (i, v) in pairs]
            @test (@inferred MA.operate!!(-, a)) === a
            @test coeffs(a) === c
            negated = [i => v for (i, v) in SA.nonzero_pairs(c)]
            @test negated == expected
            if T == BigInt
                # The `BigInt` coefficients are negated in place.
                @test all(last.(negated) .=== last.(pairs))
            end
        end
    end
    # Negation promotes `Bool` coefficients so it cannot be done in place.
    a = AlgebraElement([true, false, true], alg)
    @test MA.promote_operation(-, typeof(a)) == typeof(-a)
    b = MA.operate!!(-, a)
    @test b == -a
    @test coeffs(a) == [true, false, true]
    # Tuple storage is immutable.
    a = AlgebraElement(SA.SparseCoefficients((1, 3), (2, 3)), alg)
    b = MA.operate!!(-, a)
    @test b == -a
    @test values(coeffs(a)) == (2, 3)
end

@testset "Custom coefficient storage" begin
    alg = StarAlgebra(1.0, SA.FixedBasis{Float64,UInt32}([1.0, 2.0, 3.0, 4.0]))
    a = AlgebraElement(ACoeffs(BigInt[1, 0, 3, 0]), alg)
    b = MA.copy_if_mutable(a)
    @test typeof(b) === typeof(a)
    @test parent(b) === parent(a)
    @test SA.coeffs(b).vals !== SA.coeffs(a).vals
    MA.operate!(+, SA.coeffs(b)[1], 1)
    @test SA.coeffs(b)[1] == 2
    @test SA.coeffs(a)[1] == 1
    @test MA.operate!!(remove_leading_term, b) === b
    @test SA.coeffs(b).vals == [2, 0, 0, 0]
    @test MA.operate!!(zero, b) === b
    @test iszero(b)
    @test SA.coeffs(a).vals == [1, 0, 3, 0]

    # Allocating negation works, but missing MA support must not be hidden.
    @test (-SA.coeffs(a)).vals == [-1, 0, -3, 0]
    @test_throws MethodError MA.operate!!(-, SA.coeffs(a))
end

@testset "Key mapping" begin
    source = SA.SparseCoefficients([[2, 1], [3, 1]], BigInt[2, 3])
    output = copy(source)
    storage, coefficients = keys(output), values(output)
    shift = [1, 0]
    @test (@inferred SA.map_keys!(Base.Fix2(-, shift), output)) === output
    @test keys(output) == [[1, 1], [2, 1]]
    @test keys(output) === storage
    @test values(output) === coefficients
    @test values(output) == BigInt[2, 3]
    @test keys(source) == [[2, 1], [3, 1]]
    @test shift == [1, 0]
    empty = zero(source)
    @test (@inferred SA.map_keys!(Base.Fix2(-, shift), empty)) === empty
    @test iszero(empty)
    let output = SA.SparseCoefficients([1, 3], [2, 4])
        f = Base.Fix2(+, 1)
        SA.map_keys!(f, output)
        @test (@allocated SA.map_keys!(f, output)) == 0
        @test keys(output) == [3, 5]
    end
end

@testset "Coefficient mapping" begin
    alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0, 4.0]))
    for T in (Int, BigInt), alias in (false, true)
        for c in (
            T[2, 0, 3, 0],
            sparsevec([1, 3], T[2, 3], 4),
            SA.SparseCoefficients([1, 3], T[2, 3]),
        )
            input = AlgebraElement(c, alg)
            output = alias ? input : copy(input)
            @test @inferred(
                SA.map_coefficients_to!(output, c -> c - 2, input)
            ) === output
            @test [SA.coeffs(output)[i] for i in 1:4] == [0, 0, 1, 0]
            if !alias
                @test [SA.coeffs(input)[i] for i in 1:4] == [2, 0, 3, 0]
            end
            @test SA.map_coefficients!(c -> 2c, output; nonzero = true) ===
                  output
            @test [SA.coeffs(output)[i] for i in 1:4] == [0, 0, 2, 0]
            @test SA.map_coefficients!(zero, output) === output
            @test iszero(output)
            if c isa Union{SA.SparseCoefficients,SparseVector}
                @test isempty(collect(SA.nonzero_pairs(SA.coeffs(output))))
            end
        end
    end

    # Preserve canonical output when storage uses a different ordering.
    c = SA.SparseCoefficients([1, 2, 3], [1, 2, 3])
    out = SA.SparseCoefficients(Int[], Float64[], >)
    @test SA.map_coefficients_to!(out, c -> c / 2, c) === out
    @test keys(out) == [3, 2, 1]
    @test values(out) == [1.5, 1.0, 0.5]
    @test values(c) == [1, 2, 3]
    repeated = SA.SparseCoefficients([3, 1, 3], [1, 2, 3])
    @test SA.map_coefficients_to!(out, c -> c^2, repeated) === out
    @test keys(out) == [3, 1]
    @test values(out) == [10, 4]

    input = AlgebraElement(c, alg)
    other_alg = StarAlgebra(1.0, SA.FixedBasis([1.0, 2.0, 3.0, 5.0]))
    other = AlgebraElement(copy(c), other_alg)
    @test_throws ArgumentError SA.map_coefficients_to!(other, zero, input)
    @test parent(other) === other_alg
    @test values(SA.coeffs(other)) == [1, 2, 3]
end

end
