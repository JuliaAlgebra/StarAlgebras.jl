# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

struct DummyBasis{T} <: SA.ExplicitBasis{T,Int}
    elements::Vector{T}
end

Base.length(b::DummyBasis) = length(b.elements)
Base.getindex(b::DummyBasis, i::Int) = b.elements[i]

@testset "Basic tests" begin
    b = DummyBasis(Irrational[π, ℯ])
    a = StarAlgebra(nothing, SA.DiracMStructure(b, *))
    s(i) = sprint(show, MIME"text/plain"(), i)
    @test sprint(show, AlgebraElement([2, -1], a)) == "2·$(s(π)) - 1·$(s(ℯ))"
end
