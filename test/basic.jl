# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2021-2025: Marek Kaluba, Benoît Legat

struct PlaceholderObject end
Base.one(::PlaceholderObject) = 1

struct DummyBasis{T} <: SA.ExplicitBasis{T,Int}
    elements::Vector{T}
end

Base.length(b::DummyBasis) = length(b.elements)
Base.getindex(b::DummyBasis, i::Int) = b.elements[i]
Base.iterate(b::DummyBasis, args...) = iterate(b.elements, args...)

@testset "Basic tests" begin
    b = DummyBasis(Irrational[π, ℯ])
    a = StarAlgebra(PlaceholderObject(), b)
    s(i) = sprint(show, MIME"text/plain"(), i)
    @test sprint(show, AlgebraElement([2, -1], a)) == "2·$(s(π)) - 1·$(s(ℯ))"
    @test a == a
    @test a == StarAlgebra(PlaceholderObject(), b)
    @test a != StarAlgebra(1, b)
    b2 = DummyBasis(Irrational[π, Irrational{:γ}()])
    @test a != StarAlgebra(PlaceholderObject(), b2)

    m = SA.MappedBasis(b, float, error)
    m2 = SA.MappedBasis(b2, float, error)
    A = StarAlgebra(1.0, m)
    A2 = StarAlgebra(1.0, m2)
    @test A == StarAlgebra(1.0, m)
    @test A == StarAlgebra(1.0, SA.MappedBasis(b, float, error))
    @test A != StarAlgebra(1.0, m2)

    s = SA.SubBasis(m, Irrational[π])
    B = StarAlgebra(1.0, s)
    @test B == B
    @test B == StarAlgebra(1.0, SA.SubBasis(m, Irrational[π]))
    @test B != StarAlgebra(1.0, SA.SubBasis(m, Irrational[ℯ]))
    @test B != StarAlgebra(1.0, SA.SubBasis(m2, Irrational[π]))
end
