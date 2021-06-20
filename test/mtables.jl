@testset "TrivialMStructure" begin
    b = StarAlgebras.Basis{UInt8}(words([:a, :b, :c], radius=4))
    k = findfirst(w->length(w)==3, b)-1

    mstr = StarAlgebras.TrivialMStructure{false}(b)

    @test mstr[1, 1] == 1
    @test mstr[1, 2] == 2
    @test mstr[2, 1] == 2
    @test mstr[3, 1] == 3

    idx = b[b[2]*b[3]] # == 8
    @test mstr[2, 3] == idx

    idx = b[b[3]*b[2]]
    @test mstr[3, 2] == idx

    tmstr = StarAlgebras.TrivialMStructure{true}(b)
    @test tmstr[1, 2] == 2
    @test tmstr[2, 3] == b[StarAlgebras.star(b[2])*b[3]]
    @test tmstr[3, 2] == b[StarAlgebras.star(b[3])*b[2]]

    @test_throws StarAlgebras.ProductNotDefined mstr[k+1, k]
end

@testset "CachedMTable" begin
    b = StarAlgebras.Basis{UInt8}(words([:a, :b, :c], radius=4))
    k = findfirst(w->length(w)==3, b)-1

    mstr = StarAlgebras.CachedMTable{false}(b, table_size=(k,k))
    @test all(iszero, mstr.table)
    StarAlgebras.cache!(mstr, 1, 2)
    @test mstr.table[1, 2] == 2
    @test mstr.table[1, 1] == 0

    idx = b[b[2]*b[3]] # == 8

    @test mstr.table[2, 3] == 0
    @test mstr[2, 3] == idx
    @test mstr.table[2, 3] == idx

    @test mstr.table[1, 3] == 0
    @test mstr.table[1, 4] == 0
    StarAlgebras.cache!(mstr, [1], [3, 4])
    @test mstr.table[1, 3] == 3
    @test mstr.table[1, 4] == 4

    tmstr = StarAlgebras.CachedMTable{true}(b, table_size=(k,k))

    @test all(iszero, tmstr.table)
    @test tmstr[1, 2] == 2
    @test tmstr[2, 3] == b[StarAlgebras.star(b[2])*b[3]]
    @test tmstr[3, 2] == b[StarAlgebras.star(b[3])*b[2]]

    @test_throws StarAlgebras.ProductNotDefined mstr[k+1, k]
end

