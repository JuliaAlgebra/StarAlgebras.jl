# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

using Test
using LinearAlgebra
using SparseArrays

using StarAlgebras
import StarAlgebras as SA
import MutableArithmetics as MA

begin
    # turning GroupsCore.GroupElement into "read-only" SparseCoefficients
    import GroupsCore: GroupElement
    Base.keys(g::GroupElement) = (g,)
    Base.values(::GroupElement) = (1,)
    SA.canonical(g::GroupElement) = g
    # missing above: Base.isless for the SA.canonical
    # TODO: implement sort by hashing for non-comparable elements?

    # generalities for meaningful *-algebras
    SA.star(g::GroupElement) = inv(g)
end

for file in readdir(joinpath(@__DIR__, "..", "examples"))
    if endswith(file, ".jl")
        include(joinpath(@__DIR__, "..", "examples", file))
    end
end

# Test that `basis` implements the vector interface
# returns the same as `vector`
function vector_interface(basis, vector)
    @test eachindex(basis) == eachindex(vector)
    @test length(basis) == length(vector)
    @test size(basis) == size(vector)
    @test eltype(basis) == eltype(vector)
    @test eltype(basis) == eltype(typeof(basis))
    @test firstindex(basis) == firstindex(vector)
    @test lastindex(basis) == lastindex(vector)
    i = firstindex(basis)
    el = getindex(basis, i)
    @test getindex(basis, i) == getindex(vector, i)
end

@testset "StarAlgebras" begin
    include("basic.jl")
    # proof of concept
    using PermutationGroups
    include("perm_grp_algebra.jl")

    include("constructors.jl")
    include("group_algebra.jl")
    include("abstract_coeffs.jl")

    # free monoid algebra
    include("monoid_algebra.jl")
    include("quadratic_form.jl")

    include("caching_allocations.jl")

    # some applications:
    using Groups
    lexord(a, b) = isless(Groups.word(a), Groups.word(b))
    SA.comparable(::Type{<:Groups.FPGroupElement}) = lexord
    include("sum_of_squares.jl")
end
