# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

module TestBroadcast

using Test
import StarAlgebras as SA

struct Key
    i::Int
end
Base.isless(a::Key, b::Key) = isless(a.i, b.i)

function _test_broadcast(keys)
    c = SA.SparseCoefficients(keys, [1.0, 2.0])
    expected(values) = SA.SparseCoefficients(keys, values)
    @test 2 .* c == expected([2.0, 4.0])
    @test c .+ c == expected([2.0, 4.0])
    # Nested broadcast
    @test c .* 2 .+ 1 == expected([3.0, 5.0])
    @test 2 .* (c .+ 1) == expected([4.0, 6.0])
    # In-place broadcast
    d = similar(c)
    @test (d .= 2 .* c) === d
    @test d == expected([2.0, 4.0])
    c .= c .* c
    @test c == expected([1.0, 4.0])
    return
end

@testset "Broadcasting SparseCoefficients" begin
    # `Int` keys used to be ambiguous with `getindex(::Broadcasted, ::Integer)`
    _test_broadcast([1, 2])
    _test_broadcast([Key(1), Key(2)])
end

end
