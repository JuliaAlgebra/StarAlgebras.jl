@testset "Algebra and Elements" begin
    A = [:a, :b, :c]
    b = StarAlgebras.Basis{UInt8}(words(A, radius=2))
    l = length(b)

    RG = StarAlgebra(one(first(b)), b, (4, 4))

    a = rand(l)

    @test AlgebraElement(a, RG) isa AlgebraElement
    @test all(RG(g) isa AlgebraElement{typeof(RG)} for g in b)

    @test typeof(zero(RG)) == typeof(RG(0))
    @test typeof(one(RG)) == typeof(RG(1))

    @test_throws AssertionError AlgebraElement([1, 2, 3], RG)
    @test AlgebraElement([1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], RG) isa AlgebraElement

    x = AlgebraElement([1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], RG)
    @test AlgebraElement{Float64}(x) isa AlgebraElement
    y = AlgebraElement{Float64}(x)
    @test parent(x) === parent(y)
    @test x == y

    p = Word(A, [2, 3])
    a = RG(p)
    @test coeffs(a)[b[p]] == 1
    @test coeffs(a) isa SparseVector
    @test all(coeffs(a)[i] == 0 for i in 1:length(b) if i ≠ b[p])
    @test a(p) == 1
    @test all(a(g) == 0 for g in b if g != p)

    @test sprint(show, zero(RG)) == "0·(id)"
    @test sprint(show, one(RG)) == "1·(id)"
    @test isone(one(a))
    @test !isone(x)
    @test iszero(zero(a))
    @test sprint(show, a) == "1·b·c"
    @test sprint(show, -a) == "-1·b·c"

    @test hash(a) == hash(one(RG) + RG(p) - one(RG))

    z = zeros(l)
    z[b[p]] = 1
    @test AlgebraElement(z, RG) == a

    @test StarAlgebras.supp(a) == [p]
    @test StarAlgebras.supp_ind(a) == [b[p]]

    s = one(first(b))
    @test a(s) == 0

    a[s] = 2

    @test coeffs(a)[b[s]] == 2
    @test a[b[s]] == 2
    @test a(s) == 2

    dense_a = AlgebraElement(Vector(coeffs(a)), RG)
    @test a == dense_a
    @test hash(a) == hash(dense_a)

    @test StarAlgebras.supp_ind(a) == [b[s], b[p]] == StarAlgebras.supp_ind(dense_a)
    @test supp(a) == [s, p] == StarAlgebras.supp(dense_a)

    aa = a - RG(p)
    dense_aa = dense_a - RG(p)
    @test StarAlgebras.supp_ind(aa) == [b[s]] == StarAlgebras.supp_ind(dense_aa)
    @test supp(aa) == [s] == StarAlgebras.supp(dense_aa)

    @test sprint(show, a) == "2·(id) +1·b·c"
    @test sprint(show, -a) == "-2·(id) -1·b·c"
    z[b[s]] = 2
    @test AlgebraElement(z, RG) == a
    @test sprint(show, AlgebraElement(z, RG)) == "2.0·(id) +1.0·b·c"
    @test sprint(show, 2one(RG) - RG(p)) == "2·(id) -1·b·c"

    @test LinearAlgebra.norm(a, 1) == 3

    @test copy(a) == a
    @test copy(a) !== a
    @test coeffs(copy(a)) !== coeffs(a)
    @test parent(copy(a)) === parent(a)

    @test deepcopy(a) == a
    @test deepcopy(a) !== a
    @test coeffs(deepcopy(a)) !== coeffs(a)
    @test parent(deepcopy(a)) === parent(a)

    @testset "without basis" begin
        O = one(first(b))
        RG2 = StarAlgebra(O, RG.mstructure)
        @test_throws String zero(RG2)
        @test_throws String one(RG2)
        @test_throws String RG2(one(O))
        @test_throws String RG2(-5)

        @test sprint(show, AlgebraElement(rand(-2:2, 6), RG2)) isa String
        @test sprint(show, AlgebraElement(zeros(Int, 6), RG2)) isa String
    end
end
