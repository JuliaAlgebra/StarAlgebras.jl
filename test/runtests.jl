using StarAlgebras
using Test
using LinearAlgebra
using SparseArrays

include("test_example_words.jl")

@testset "StarAlgebras" begin
    include("mtables.jl")
    include("constructors.jl")

    using GroupsCore
    StarAlgebras.star(g::GroupsCore.GroupElement) = inv(g)

    using SparseArrays
    if VERSION < v"1.9"
        Base.sum(v::SparseVector) = sum(nonzeros(v))
    end

    using PermutationGroups
    include("arithmetic.jl")

    using Groups
    include("sum_of_squares.jl")
end
