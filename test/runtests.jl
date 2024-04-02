using Test
using Random
using LinearAlgebra
using SparseArrays

using StarAlgebras
import StarAlgebras as SA
import MutableArithmetics as MA

begin
    # turning GroupsCore.GroupElement into SparseCoefficients
    import GroupsCore: GroupElement
    Base.iszero(::GroupElement) = false
    Base.keys(g::GroupElement) = (g,)
    Base.values(::GroupElement) = (1,)
    SA.canonical(g::GroupElement) = g
    function Base.getindex(g::GroupElement, h::GroupElement)
        return ifelse(g == h, SA.aug(g), 0 * SA.aug(g))
    end
    # missing above: Base.isless for the SA.canonical
    # TODO: implement sort by hashing for non-comparable elements?

    # generalities for meaningful *-algebras
    SA.star(g::GroupElement) = inv(g)
    SA.aug(::GroupElement) = 1
end

# implementations of the free monoid over an alphabet
include("test_example_words.jl")

@testset "StarAlgebras" begin
    # proof of concept
    using PermutationGroups
    include("perm_grp_algebra.jl")

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

