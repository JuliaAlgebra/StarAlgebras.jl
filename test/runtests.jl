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
end

# implementations of the free monoid over an alphabet
include("test_example_words.jl")

@testset "StarAlgebras" begin
    # proof of concept
    using PermutationGroups
    include("perm_grp_algebra.jl")

    # arithmetic for perm group algebra and the free monoid algebra
    include("arithmetic.jl")

    # include("mtables.jl")
    include("constructors.jl")

    # using SparseArrays
    # if VERSION < v"1.9"
    #     Base.sum(v::SparseVector) = sum(nonzeros(v))
    # end

    # some applications:
    using Groups
    function Base.isless(g::Groups.FPGroupElement, h::Groups.FPGroupElement)
        return isless(Groups.word(g), Groups.word(h))
    end
    include("sum_of_squares.jl")
end

