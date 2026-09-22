# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

# Reserve destination capacity at the call site. Restore inputs after warmup
# so aliasing and cancellation still exercise the full merge when measured.
function test_allocation_free_merge!(
    run!::F,
    output,
    expected,
    storage,
) where {F}
    original = map(copy, storage)
    @test @inferred(run!()) === output
    @test output == expected
    for (v, saved) in zip(storage, original)
        resize!(v, length(saved))
        copyto!(v, saved)
    end
    @test @allocated(run!()) == 0
    @test output == expected
    return
end

function test_coefficient_mapping_allocations(f::F, alias, nonzero, n) where {F}
    input = SA.SparseCoefficients(collect(1:n), collect(1:n))
    original = copy(input)
    output = alias ? input : zero(input)
    sizehint!(keys(output), n)
    sizehint!(values(output), n)
    @test @inferred(SA.map_coefficients_to!(output, f, input; nonzero)) ===
          output
    for (v, saved) in
        ((keys(input), keys(original)), (values(input), values(original)))
        resize!(v, n)
        copyto!(v, saved)
    end
    @test @allocated(SA.map_coefficients_to!(output, f, input; nonzero)) == 0
    @test collect(SA.nonzero_pairs(output)) ==
          [i => f(i) for i in 1:n if !iszero(f(i))]
    return
end

@testset "Allocation-free coefficient mapping" begin
    for n in (0, 100), alias in (false, true)
        test_coefficient_mapping_allocations(Base.Fix2(mod, 3), alias, false, n)
        test_coefficient_mapping_allocations(-, alias, true, n)
    end
end

@testset "Allocation-free sorted merge" begin
    for n in (0, 1, 100),
        rev in (false, true),
        alias in (:none, :left, :right, :both),
        cancel in (false, true)

        x = collect(1:2:2n)
        y = alias === :both ? x : cancel ? copy(x) : collect(2:3:2n)
        sort!(x; rev)
        sort!(y; rev)
        output =
            alias === :left || alias === :both ? x :
            alias === :right ? y : Int[]
        sizehint!(output, length(x) + length(y))
        combine = cancel ? (-) : SA.first_of
        expected = cancel ? Int[] : sort!(unique([x; y]); rev)
        test_allocation_free_merge!(
            () -> SA.merge_sorted!(
                output,
                x,
                y;
                lt = isless,
                combine,
                filter = !iszero,
                rev,
            ),
            output,
            expected,
            (x, y),
        )
    end
end

@testset "Allocation-free coefficient merge" begin
    for T in (Int, Float64),
        n in (0, 1, 100),
        alias in (:none, :left, :right, :both)

        x = SA.SparseCoefficients(collect(0:(n-1)), ones(T, n))
        y =
            alias === :both ? x :
            SA.SparseCoefficients(collect(1:n), -ones(T, n))
        expected = x + y
        output =
            alias === :left || alias === :both ? x :
            alias === :right ? y : zero(x)
        sizehint!(keys(output), 2n)
        sizehint!(values(output), 2n)
        storage = (keys(x), values(x), keys(y), values(y))
        # This view must expose its pair type to Reverse without a Generator.
        pairs = Iterators.Reverse(SA._CoefficientPairs(x))
        @test eltype(pairs) === Pair{Int,T}
        @test Base.IteratorEltype(typeof(pairs)) === Base.HasEltype()
        test_allocation_free_merge!(
            () -> MA.operate_to!(output, +, x, y),
            output,
            expected,
            storage,
        )
    end
    for T in (Int, Float64)
        x = SA.SparseCoefficients(collect(1:100), ones(T, 100))
        y = copy(x)
        expected = x + y
        sizehint!(keys(x), 200)
        sizehint!(values(x), 200)
        test_allocation_free_merge!(
            () -> MA.operate!!(+, x, y),
            x,
            expected,
            (keys(x), values(x), keys(y), values(y)),
        )
    end
end

@testset "Allocation-free ordered term accumulation" begin
    # Tuple keys and machine coefficients have nonallocating scalar products.
    for T in (Int, Float64),
        n in (0, 1, 100),
        op in (MA.add_mul, MA.sub_mul),
        left in (true, false),
        alias in (false, true),
        shift in (0, 1)

        c = SA.SparseCoefficients([(i, 0) for i in 0:(n-1)], ones(T, n), grlex)
        f = AlgebraElement(c, biv_alg)
        g = alias ? f : AlgebraElement(copy(c), biv_alg)
        t = biv_term(one(T), (shift, 0))
        a, b = left ? (t, g) : (g, t)
        expected = op(f, a, b)
        sizehint!(keys(c), 2n)
        sizehint!(values(c), 2n)
        test_allocation_free_merge!(
            () -> MA.operate!!(op, f, a, b),
            f,
            expected,
            (keys(c), values(c), keys(SA.coeffs(g)), values(SA.coeffs(g))),
        )
    end
end
