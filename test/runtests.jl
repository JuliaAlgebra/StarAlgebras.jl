using StarAlgebras
using Test
using LinearAlgebra
using SparseArrays

include("test_example_words.jl")

@testset "StarAlgebras" begin
    include("mtables.jl")
    include("constructors.jl")
    include("arithmetic.jl")

    if VERSION >= v"1.3.0"
        using Pkg
        Pkg.add("Groups")
        include("sum_of_squares.jl")
    end
end
