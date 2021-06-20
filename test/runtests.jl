using StarAlgebras
using Test
using LinearAlgebra
using SparseArrays

include("test_example_words.jl")

@testset "StarAlgebras" begin
   include("mtables.jl")
   include("constructors.jl")
   include("arithmetic.jl")
end
