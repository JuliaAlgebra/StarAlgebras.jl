using Test, LinearAlgebra, SparseArrays, PermutationGroups, StarAlgebras
import StarAlgebras as SA
import MutableArithmetics as MA
import GroupsCore: GroupElement

# Treat group elements as read-only coefficient containers.
Base.keys(g::GroupElement) = (g,)
Base.values(::GroupElement) = (1,)
SA.canonical(g::GroupElement) = g
SA.star(g::GroupElement) = inv(g)

if !@isdefined aug
    include("../examples/diracs_augmented.jl")
end
