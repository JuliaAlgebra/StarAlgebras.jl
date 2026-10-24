# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

using Test
using LinearAlgebra
using SparseArrays

using StarAlgebras
import StarAlgebras as SA
import MutableArithmetics as MA

for file in readdir(joinpath(@__DIR__, "..", "examples"))
    if endswith(file, ".jl")
        include(joinpath(@__DIR__, "..", "examples", file))
    end
end

include("group_utils.jl")

# Test that `basis` implements the vector interface
# returns the same as `vector`
function test_vector_interface(basis, vector = collect(basis))
    @test length(basis) == length(vector)
    @test size(basis) == size(vector)
    @test eltype(basis) == eltype(vector)
    @test eltype(basis) == eltype(typeof(basis))
    @test last(basis) == last(vector)
    i = lastindex(basis)
    @test i == lastindex(vector)
    @test getindex(basis, i) == last(vector)
    @test first(basis) == first(vector)
    i = firstindex(basis)
    @test i == firstindex(vector)
    @test getindex(basis, i) == first(vector)
end

@testset "StarAlgebras" begin
    include("term.jl")
    include("remove_leading_term.jl")
    include("mutable_arithmetics.jl")
    include("term_arithmetic.jl")
    include("term_addition.jl")
    include("scalar_arithmetic.jl")
    include("extensions.jl")
    include("parent_promotion.jl")
    include("zero_for.jl")
    include("basis_conversion.jl")
    include("merge_allocations.jl")
    include("basic.jl")
    include("broadcast.jl")
    # proof of concept
    using PermutationGroups
    include("merge_sorted.jl")
    include("perm_grp_algebra.jl")

    include("constructors.jl")
    include("group_algebra.jl")
    include("abstract_coeffs.jl")

    # free monoid algebra
    include("monoid_algebra.jl")
    include("quadratic_form.jl")

    include("caching_allocations.jl")
    include("graded_lex.jl")

    # some applications:
    using Groups
    lexord(a, b) = isless(Groups.word(a), Groups.word(b))
    SA.comparable(::Type{<:Groups.FPGroupElement}) = lexord
    include("sum_of_squares.jl")

    include("invalidations.jl")
end
