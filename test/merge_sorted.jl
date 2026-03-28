@testset "merge_sorted" begin
    dedup_kw = (; combine = SA.first_of, filter = _ -> true)
    @testset "Vector dedup merge" begin
        # Disjoint sorted vectors
        @test SA.merge_sorted([1, 3, 5], [2, 4, 6]; lt = isless, dedup_kw...) ==
              [1, 2, 3, 4, 5, 6]

        # Overlapping sorted vectors
        @test SA.merge_sorted([1, 3, 5], [3, 5, 7]; lt = isless, dedup_kw...) ==
              [1, 3, 5, 7]

        # Identical vectors
        @test SA.merge_sorted([1, 2, 3], [1, 2, 3]; lt = isless, dedup_kw...) ==
              [1, 2, 3]

        # One empty
        @test SA.merge_sorted(Int[], [1, 2, 3]; lt = isless, dedup_kw...) ==
              [1, 2, 3]
        @test SA.merge_sorted([1, 2, 3], Int[]; lt = isless, dedup_kw...) ==
              [1, 2, 3]

        # Both empty
        @test SA.merge_sorted(Int[], Int[]; lt = isless, dedup_kw...) == Int[]

        # Single elements
        @test SA.merge_sorted([1], [2]; lt = isless, dedup_kw...) == [1, 2]
        @test SA.merge_sorted([2], [1]; lt = isless, dedup_kw...) == [1, 2]
        @test SA.merge_sorted([1], [1]; lt = isless, dedup_kw...) == [1]

        # Custom lt (reverse order)
        @test SA.merge_sorted([5, 3, 1], [4, 2]; lt = >, dedup_kw...) ==
              [5, 4, 3, 2, 1]
        @test SA.merge_sorted([5, 3, 1], [5, 2]; lt = >, dedup_kw...) ==
              [5, 3, 2, 1]

        # Mixed element types
        @test SA.merge_sorted([1, 3], [2.0, 4.0]; lt = isless, dedup_kw...) ==
              [1.0, 2.0, 3.0, 4.0]
    end

    @testset "Vector combine merge" begin
        # Combine duplicate keys by summing
        v1 = [1, 3, 5]
        v2 = [3, 5, 7]
        result = SA.merge_sorted(
            v1,
            v2;
            lt = isless,
            combine = +,
            filter = x -> !iszero(x),
        )
        @test result == [1, 6, 10, 7]

        # Filter zeros
        v1 = [1, 2, 3]
        v2 = [2, 3, 4]
        result = SA.merge_sorted(
            v1,
            v2;
            lt = isless,
            combine = +,
            filter = x -> x < 6,
        )
        @test result == [1, 4, 4]  # 2+2=4 kept, 3+3=6 filtered, 4 kept

        # No filter (keep all)
        result = SA.merge_sorted(
            v1,
            v2;
            lt = isless,
            combine = +,
            filter = _ -> true,
        )
        @test result == [1, 4, 6, 4]
    end

    @testset "In-place merge_sorted!" begin
        result = Vector{Int}(undef, 6)
        SA.merge_sorted!(
            result,
            [1, 3, 5],
            [2, 4, 6];
            lt = isless,
            combine = +,
            filter = x -> !iszero(x),
        )
        @test result == [1, 2, 3, 4, 5, 6]

        result = Vector{Int}(undef, 6)
        SA.merge_sorted!(
            result,
            [1, 3, 5],
            [3, 5, 7];
            lt = isless,
            combine = +,
            filter = x -> !iszero(x),
        )
        @test result == [1, 6, 10, 7]
    end

    @testset "Tuple dedup merge" begin
        # Disjoint
        @test SA.merge_sorted((1, 3, 5), (2, 4, 6); lt = isless, dedup_kw...) ==
              (1, 2, 3, 4, 5, 6)

        # Overlapping
        @test SA.merge_sorted((1, 3, 5), (3, 5, 7); lt = isless, dedup_kw...) ==
              (1, 3, 5, 7)

        # Identical
        @test SA.merge_sorted((1, 2, 3), (1, 2, 3); lt = isless, dedup_kw...) ==
              (1, 2, 3)

        # One empty
        @test SA.merge_sorted((), (1, 2, 3); lt = isless, dedup_kw...) ==
              (1, 2, 3)
        @test SA.merge_sorted((1, 2, 3), (); lt = isless, dedup_kw...) ==
              (1, 2, 3)

        # Both empty
        @test SA.merge_sorted((), (); lt = isless, dedup_kw...) == ()

        # Single elements
        @test SA.merge_sorted((1,), (2,); lt = isless, dedup_kw...) == (1, 2)
        @test SA.merge_sorted((2,), (1,); lt = isless, dedup_kw...) == (1, 2)
        @test SA.merge_sorted((1,), (1,); lt = isless, dedup_kw...) == (1,)

        # Custom lt (reverse order)
        @test SA.merge_sorted((5, 3, 1), (4, 2); lt = >, dedup_kw...) ==
              (5, 4, 3, 2, 1)
        @test SA.merge_sorted((5, 3, 1), (5, 2); lt = >, dedup_kw...) ==
              (5, 3, 2, 1)
    end

    @testset "Vector rev keyword" begin
        # rev=true with default isless reverses the order
        @test SA.merge_sorted([5, 3, 1], [4, 2]; lt = isless, rev = true, dedup_kw...) ==
              [5, 4, 3, 2, 1]
        @test SA.merge_sorted([5, 3, 1], [5, 2]; lt = isless, rev = true, dedup_kw...) ==
              [5, 3, 2, 1]

        # rev=true with combine
        result = SA.merge_sorted([5, 3, 1], [3, 1]; lt = isless, rev = true, combine = +, filter = _ -> true)
        @test result == [5, 6, 2]

        # rev=false is the default behavior
        @test SA.merge_sorted([1, 3, 5], [2, 4, 6]; lt = isless, rev = false, dedup_kw...) ==
              [1, 2, 3, 4, 5, 6]
    end

    @testset "In-place rev keyword" begin
        result = Vector{Int}(undef, 5)
        SA.merge_sorted!(result, [5, 3, 1], [4, 2]; lt = isless, rev = true, combine = SA.first_of, filter = _ -> true)
        @test result == [5, 4, 3, 2, 1]
    end

    @testset "Tuple combine merge" begin
        @test SA.merge_sorted(
            (1, 3, 5),
            (3, 5, 7);
            lt = isless,
            combine = +,
            filter = _ -> true,
        ) == (1, 6, 10, 7)
        # Filter
        @test SA.merge_sorted(
            (1, 3, 5),
            (3, 5, 7);
            lt = isless,
            combine = +,
            filter = x -> x < 10,
        ) == (1, 6, 7)
    end

    @testset "Tuple rev keyword" begin
        @test SA.merge_sorted((5, 3, 1), (4, 2); lt = isless, rev = true, dedup_kw...) ==
              (5, 4, 3, 2, 1)
        @test SA.merge_sorted((5, 3, 1), (5, 2); lt = isless, rev = true, dedup_kw...) ==
              (5, 3, 2, 1)

        # rev=true with combine
        @test SA.merge_sorted((5, 3, 1), (3, 1); lt = isless, rev = true, combine = +, filter = _ -> true) ==
              (5, 6, 2)
    end

    @testset "multi_findsorted" begin
        # All found
        @test SA.multi_findsorted([1, 2, 3], [1, 2, 3]) == [1, 2, 3]

        # Superset x, subset y
        @test SA.multi_findsorted([1, 2, 3, 4, 5], [2, 4]) == [0, 1, 0, 2, 0]

        # Subset x, superset y
        @test SA.multi_findsorted([2, 4], [1, 2, 3, 4, 5]) == [2, 4]

        # Disjoint
        @test SA.multi_findsorted([1, 3, 5], [2, 4, 6]) == [0, 0, 0]

        # Empty x
        @test SA.multi_findsorted(Int[], [1, 2, 3]) == Int[]

        # Empty y
        @test SA.multi_findsorted([1, 2, 3], Int[]) == [0, 0, 0]

        # Custom lt (reverse order)
        @test SA.multi_findsorted([5, 3, 1], [5, 4, 3, 2, 1]; lt = >) ==
              [1, 3, 5]
    end

    @testset "merge_bases_with_maps" begin
        using PermutationGroups
        G = PermGroup(perm"(1,2,3)", perm"(1,2)")
        db = SA.DiracBasis(G)
        elts = sort(collect(G))

        b1 = SA.SubBasis(db, elts[1:3])
        b2 = SA.SubBasis(db, elts[2:5])
        @test b1.is_sorted
        @test b2.is_sorted

        merged, I1, I2 = SA.merge_bases_with_maps(b1, b2)
        @test length(merged) == 5
        @test merged.keys == elts[1:5]

        # Check index mapping: I1[i] = index of merged[i] in b1, 0 if absent
        for i in eachindex(I1)
            if I1[i] != 0
                @test merged.keys[i] == b1.keys[I1[i]]
            end
        end
        for i in eachindex(I2)
            if I2[i] != 0
                @test merged.keys[i] == b2.keys[I2[i]]
            end
        end

        # Disjoint bases
        b3 = SA.SubBasis(db, elts[1:2])
        b4 = SA.SubBasis(db, elts[4:6])
        merged2, I3, I4 = SA.merge_bases_with_maps(b3, b4)
        @test length(merged2) == 5
        @test all(==(0), I3[3:5])
        @test all(==(0), I4[1:2])
        @test merged2 == SA.merge_bases(b3, b4)

        # Identical bases
        merged3, I5, I6 = SA.merge_bases_with_maps(b1, b1)
        @test merged3.keys == b1.keys
        @test I5 == I6 == collect(1:length(b1))
    end
end
