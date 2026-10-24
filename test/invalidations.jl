# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

using Test

# Invalidations can only be measured when the package is loaded for the first
# time so we run the measurement in a fresh process. The dependencies (and the
# triggers of the extensions) are loaded first so that only the invalidations
# caused by methods defined in StarAlgebras (or its extensions) are counted.
const INVALIDATIONS_SCRIPT = """
using LinearAlgebra, SparseArrays, MutableArithmetics, MathOptInterface, JuMP
using SnoopCompileCore
invs = @snoop_invalidations using StarAlgebras
using SnoopCompile
dir = pkgdir(StarAlgebras)
trees = filter(invalidation_trees(invs)) do tree
    return startswith(String(tree.method.file), dir)
end
show(stderr, MIME"text/plain"(), trees)
println(stderr)
println(sum(SnoopCompile.countchildren, trees; init = 0))
"""

@testset "Invalidations" begin
    project = dirname(Base.active_project())
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$project -e $INVALIDATIONS_SCRIPT`
    count = parse(Int, readchomp(cmd))
    @test count == 0
end
