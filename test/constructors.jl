struct CustomLaTeXPrint
    s::String
end

Base.:-(s::CustomLaTeXPrint) = s
Base.iszero(::CustomLaTeXPrint) = false
function Base.show(io::IO, ::MIME"text/latex", s::CustomLaTeXPrint)
    return print(io, s.s)
end

@testset "Algebra and Elements" begin
    alph = [:a, :b, :c]
    A★ = FreeWords(alph)
    B = SA.DiracBasis(A★)
    RG = StarAlgebra(A★, B)
    @test typeof(@inferred basis(RG)) == MA.promote_operation(basis, typeof(RG))

    @test typeof(zero(RG)) == typeof(RG(0))
    @test typeof(one(RG)) == typeof(RG(1))

    @test one(RG) == RG(1)
    @test RG(1) == one(RG(0))

    @test zero(RG) == RG(0)
    @test RG(0) == zero(RG(1))

    x = let A = alph, n = 7
        elts = [Word(A, rand(1:length(A), rand(0:7))) for _ = 1:n]
        vals = rand(-5:5, n)
        coeffs = SA.SparseCoefficients(elts, vals)
        MA.operate!(SA.canonical, coeffs)
    end

    @test AlgebraElement(x, RG) isa AlgebraElement

    X = AlgebraElement(x, RG)
    @test typeof(@inferred coeffs(X)) == MA.promote_operation(coeffs, typeof(X))
    @test typeof(@inferred basis(X)) == MA.promote_operation(basis, typeof(X))

    @test AlgebraElement{Float64}(X) isa AlgebraElement
    Y = AlgebraElement{Float64}(X)
    @test parent(X) === parent(Y)
    @test X == Y

    p = Word(alph, [2, 3])
    a = RG(p)
    @test coeffs(a)[B[p]] == 1
    @test coeffs(a) isa SA.SparseCoefficients
    @test length(collect(SA.nonzero_pairs(coeffs(a)))) == 1
    @test a(p) == 1

    let b = collect(Iterators.take(B, 121))
        @test all(coeffs(a)[x] == 0 for x in b if x ≠ p)
        @test all(a(g) == 0 for g in b if g != p)
    end

    @test sprint(show, zero(RG)) == "0·(id)"
    @test sprint(show, one(RG)) == "1·(id)"
    @test isone(one(a))
    @test !isone(a)
    @test iszero(zero(a))
    @test sprint(show, a) == "1·b·c"
    @test sprint(show, -a) == "-1·b·c"

    @test hash(a) == hash(one(RG) + RG(p) - one(RG))

    @test SA.supp(a) == [p]
    @test collect(SA.nonzero_pairs(coeffs(a))) == [(p => 1)]

    s = one(first(B))
    @test a(s) == 0

    a[s] = 2

    @test coeffs(a)[s] == 2
    @test a(s) == 2

    # dense_a = AlgebraElement(Vector(coeffs(a)), RG)
    # @test a == dense_a
    # @test hash(a) == hash(dense_a)

    # @test SA.supp_ind(a) == [b[s], b[p]] == SA.supp_ind(dense_a)
    # @test SA.supp(a) == [s, p] == SA.supp(dense_a)

    aa = a - RG(p)
    # dense_aa = dense_a - RG(p)
    # @test SA.supp_ind(aa) == [b[s]] == SA.supp_ind(dense_aa)
    @test SA.supp(aa) == [s]

    @test sprint(show, a) == "2·(id) + 1·b·c"
    @test sprint(show, -a) == "-2·(id) - 1·b·c"
    Z = AlgebraElement{Float64}(a)
    @test Z == a
    @test sprint(show, Z) == "2.0·(id) + 1.0·b·c"
    @test sprint(show, 2one(RG) - RG(p)) == "2·(id) - 1·b·c"
    @test sprint(show, (2 + im) * one(RG) - (3im) * RG(p)) ==
          "(2 + 1im)·(id) + (0 - 3im)·b·c"

    @test sprint(print, (2 + im) * one(RG) - (3im) * RG(p)) ==
          "(2 + 1im)·(id) + (0 - 3im)·b·c"
    @test sprint(show, 1e-9 * one(RG)) == "1.0e-9·(id)"
    @test sprint((io, x) -> show(io, "text/latex", x), 1e-9 * one(RG)) ==
          "\$\$ 1.0 \\cdot 10^{-9} \\cdot (id) \$\$"

    @test LinearAlgebra.norm(a, 1) == 3

    @test copy(a) == a
    @test copy(a) !== a
    @test coeffs(copy(a)) !== coeffs(a)
    @test parent(copy(a)) === parent(a)

    @test deepcopy(a) == a
    @test deepcopy(a) !== a
    @test coeffs(deepcopy(a)) !== coeffs(a)
    @test parent(deepcopy(a)) === parent(a)

    latex = CustomLaTeXPrint(" \$\$ \\[\\(α_β∀ \\) \\]\t  \$\$")
    a = SA.AlgebraElement(SA.SparseCoefficients([p], [latex]), RG)
    # Tests that `getindex` works even if `zero(typeof(latex))` is not defined
    @test SA.coeffs(a)[p] == latex
    @test sprint((io, x) -> show(io, "text/latex", x), a) == "\$\$ (α_β∀) \\cdot b·c \$\$"
    # Test that the check for `\\)` handles unicode well
    latex = CustomLaTeXPrint("\\(β∀")
    @test sprint(
        (io, x) -> show(io, "text/latex", x),
        SA.AlgebraElement(SA.SparseCoefficients([p], [latex]), RG),
    ) == "\$\$ (\\(β∀) \\cdot b·c \$\$"

end
