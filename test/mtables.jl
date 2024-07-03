@testset "TrivialMStructure" begin
    b = StarAlgebras.Basis{UInt8}(words([:a, :b, :c]; radius = 4))

    mstr = StarAlgebras.TrivialMStructure(b)

    @test mstr[1, 1] == 1
    @test mstr[1, 2] == 2
    @test mstr[-1, 2] == 2
    @test mstr[2, 1] == 2
    @test mstr[-2, 1] == 3
    @test mstr[3, 1] == 3
    @test mstr[-3, 1] == 2

    idx = b[b[2]*b[3]] # == 8
    @test mstr[2, 3] == idx

    idx = b[b[3]*b[2]]
    @test mstr[3, 2] == idx
    @test mstr[-2, 3] == b[star(b[2])*b[3]]
    @test mstr[-3, 2] == b[star(b[3])*b[2]]

    @test sprint(show, MIME"text/plain"(), mstr) ==
          "TrivialMStructure over basis with $(length(basis(mstr))) elements"

    k = findfirst(w -> length(w) == 3, b)
    @test_throws StarAlgebras.ProductNotDefined mstr[k, k-1]
    @test_throws StarAlgebras.ProductNotDefined mstr[k-1, k]

    try
        w = mstr[k, k-1]
    catch ex
        @test ex isa StarAlgebras.ProductNotDefined
        @test sprint(Base.showerror, ex) ==
              "Product of elements 14 and 13 is not defined on the basis or the multiplicative structure could not be completed: a·a·a · c·c = a·a·a·c·c."
    end
end

@testset "MTable" begin
    b = StarAlgebras.Basis{UInt16}(words([:a, :b, :c, :d]; radius = 4))
    k = findfirst(w -> length(w) == 3, b) - 1
    mstr = StarAlgebras.MTable(b; size = (k, k))

    @test_throws String StarAlgebras.basis(mstr)

    @test mstr isa StarAlgebras.MTable{UInt16}
    @test all(mstr[i, j] == b[b[i]*b[j]] for i in axes(mstr, 1) for j in axes(mstr, 2))
    @test all(
        mstr[-i, j] == b[star(b[i])*b[j]] for i in axes(mstr, 1) for j in axes(mstr, 2)
    )
    @test all(
        mstr[i, -j] == b[b[i]*star(b[j])] for i in axes(mstr, 1) for j in axes(mstr, 2)
    )
    @test all(
        mstr[-i, -j] == b[star(b[i])*star(b[j])] for i in axes(mstr, 1) for
        j in axes(mstr, 2)
    )
end

@testset "CachedMTable" begin
    b = StarAlgebras.Basis{UInt8}(words([:a, :b, :c]; radius = 4))
    k = findfirst(w -> length(w) == 3, b) - 1

    mstr = StarAlgebras.CachedMTable(b; table_size = (k, k))
    @test mstr isa StarAlgebras.CachedMTable{UInt8,Word{Symbol}}
    @test mstr.table.table isa Matrix{UInt8}

    @test_throws StarAlgebras.ProductNotDefined StarAlgebras._check(
        mstr.table.table,
        StarAlgebras.basis(mstr),
    )

    StarAlgebras.complete!(mstr)
    @test all(!iszero, mstr.table)

    mstr_sparse = StarAlgebras.CachedMTable(b, spzeros(UInt8, k, k))
    @test issparse(mstr_sparse.table.table)

    StarAlgebras.complete!(mstr_sparse.table.table, basis(mstr))
    @test all(!iszero, mstr.table.table)

    @test mstr == mstr_sparse

    mstr = StarAlgebras.CachedMTable(b; table_size = (k, k))

    @test all(iszero, mstr.table.table)
    StarAlgebras.cache!(mstr, 1, 2)
    @test mstr.table.table[1, 2] == 2
    @test mstr.table.table[1, 1] == 0

    idx = b[b[2]*b[3]]
    @test mstr.table.table[2, 3] == 0
    @test_throws StarAlgebras.ProductNotDefined mstr.table[2, 3]
    @test mstr[2, 3] == idx
    @test mstr.table[2, 3] == idx
    @test mstr.table.table[2, 3] == idx

    @test mstr.table.table[1, 3] == 0
    @test mstr.table.table[1, 4] == 0
    StarAlgebras.cache!(mstr, [1], [3, 4])
    @test mstr.table[1, 3] == 3
    @test mstr.table[1, 4] == 4

    @test_throws StarAlgebras.ProductNotDefined mstr[k+1, k]

    mstr = StarAlgebras.CachedMTable(b; table_size = (k, k))
    @test all(iszero, mstr.table.table)
    @test mstr[-1, 2] == 2
    @test mstr[-2, 3] == b[star(b[2])*b[3]]
    @test mstr[3, -2] == b[b[3]*star(b[2])]
    @test mstr[-3, -2] == b[star(b[3])*star(b[2])]
end
