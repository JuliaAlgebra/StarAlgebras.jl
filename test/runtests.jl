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

# implementations of the free monoid over an alphabet
include("test_example_words.jl")

@testset "StarAlgebras" begin
    # proof of concept
    using PermutationGroups
    include("perm_grp_algebra.jl")

    include("constructors.jl")
    include("group_algebra.jl")
    include("abstract_coeffs.jl")

    # free monoid algebra
    include("monoid_algebra.jl")

    include("caching_allocations.jl")

    # some applications:
    using Groups
    function Base.isless(g::Groups.FPGroupElement, h::Groups.FPGroupElement)
        return isless(Groups.word(g), Groups.word(h))
    end
    include("sum_of_squares.jl")
end

