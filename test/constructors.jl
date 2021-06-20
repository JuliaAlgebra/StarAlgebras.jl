@testset "Algebra and Elements Constructors" begin
    A = [:a, :b, :c]
    b = StarAlgebras.Basis{UInt8}(words(A, radius=2))
    l = length(b)

    RG = StarAlgebras.StarAlgebra(Symbol, b, (4, 4))

    a = rand(l)

    @test StarAlgebras.AlgebraElement(a, RG) isa StarAlgebras.AlgebraElement
    @test all(RG(g) isa StarAlgebras.AlgebraElement{typeof(RG)} for g in b)

    @test_throws AssertionError StarAlgebras.AlgebraElement([1,2,3], RG)
    @test StarAlgebras.AlgebraElement([1,2,3,0,0,0,0,0,0,0,0,0,0], RG) isa StarAlgebras.AlgebraElement

    p = Word(A, [2,3])
    a = RG(p)
    @test StarAlgebras.coeffs(a) isa SparseVector
    @test StarAlgebras.coeffs(a)[b[p]] == 1
    @test all(StarAlgebras.coeffs(a)[i] == 0 for i in 1:length(b) if i ≠ b[p])
    @test a(p) == 1
    @test all(a(g) == 0 for g in b if g != p)

    @test sprint(show, a) == "1·b·c"
    @test sprint(show, -a) == "-1·b·c"

    z = zeros(l)
    z[b[p]] = 1
    @test StarAlgebras.AlgebraElement(z, RG) == a

    @test StarAlgebras.supp(a) == [p]
    @test StarAlgebras.supp_ind(a) == [b[p]]

    s = one(first(b))
    @test a(s) == 0

    a[s] = 2

    @test StarAlgebras.coeffs(a)[b[s]] == 2
    @test a[b[s]] == 2
    @test a(s) == 2

    @test StarAlgebras.supp(a) == [s, p]
    @test StarAlgebras.supp_ind(a) == [b[s], b[p]]

    @test sprint(show, a) == "2·(id) +1·b·c"
    @test sprint(show, -a) == "-2·(id) -1·b·c"
    z[b[s]] = 2
    @test StarAlgebras.AlgebraElement(z, RG) == a
    @test sprint(show, StarAlgebras.AlgebraElement(z, RG)) == "2.0·(id) +1.0·b·c"
end
