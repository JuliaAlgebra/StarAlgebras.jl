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
    @test tmstr[2, 3] == b[star(b[2])*b[3]]
    @test tmstr[3, 2] == b[star(b[3])*b[2]]

    @test_throws StarAlgebras.ProductNotDefined mstr[k+1, k]

    try
        w = mstr[k+1, k]
    catch ex
        @test ex isa StarAlgebras.ProductNotDefined
        @test sprint(Base.showerror, ex) == "Product of elements 14 and 13 is not defined on the basis or the multiplicative structure could not be completed: a·a·a · c·c = a·a·a·c·c."
    end
end

@testset "MTable" begin
    b = StarAlgebras.Basis{UInt16}(words([:a, :b, :c, :d], radius=4))
    k = findfirst(w->length(w)==3, b)-1

    mstr = StarAlgebras.MTable(b, table_size=(k, k))

    @test mstr isa StarAlgebras.MTable{UInt16, false}

    @test all(mstr[i,i]≠1 for i in 2:size(mstr, 1))
    @test all(mstr[1,i]==i for i in 1:size(mstr, 2))
    @test all(mstr[i,1]==i for i in 1:size(mstr, 1))

    tmstr = StarAlgebras.MTable{true}(b, table_size=(k, k))

    @test tmstr isa StarAlgebras.MTable{UInt16, true}
    @test all(tmstr[i,i]!=1 for i in 2:size(tmstr, 1))
    @test all(tmstr[1,i]==i for i in 1:size(tmstr, 2))
    @test all(tmstr[i,1]≠ i for i in 1:size(tmstr, 1) if b[i] != star(b[i]))
end

@testset "CachedMTable" begin
    b = StarAlgebras.Basis{UInt8}(words([:a, :b, :c], radius = 4))
    k = findfirst(w -> length(w) == 3, b) - 1

    @test StarAlgebras.CachedMTable(b, table_size = (k, k)) isa
          StarAlgebras.CachedMTable{Word{Symbol},UInt8,typeof(b),Matrix{UInt8},false}

    @test StarAlgebras.CachedMTable{true}(b, table_size = (k, k)) isa
          StarAlgebras.CachedMTable{Word{Symbol},UInt8,typeof(b),Matrix{UInt8},true}

    @test StarAlgebras.CachedMTable(b, spzeros(UInt8, k, k)) isa StarAlgebras.CachedMTable{
        Word{Symbol},
        UInt8,
        typeof(b),
        SparseMatrixCSC{UInt8,Int64},
        false,
    }

    @test StarAlgebras.CachedMTable{true}(b, spzeros(UInt8, k, k)) isa
          StarAlgebras.CachedMTable{
        Word{Symbol},
        UInt8,
        typeof(b),
        SparseMatrixCSC{UInt8,Int64},
        true,
    }

    for mstr in [
        StarAlgebras.CachedMTable{false}(b, table_size = (k, k)),
        StarAlgebras.CachedMTable{true}(b, spzeros(UInt8, k, k)),
    ]

        @test all(iszero, mstr.table)
        StarAlgebras.cache!(mstr, 1, 2)
        @test mstr.table[1, 2] == 2
        @test mstr.table[1, 1] == 0

        idx = StarAlgebras._istwisted(mstr) ? b[star(b[2])*b[3]] : b[b[2]*b[3]]

        @test mstr.table[2, 3] == 0
        @test mstr[2, 3] == idx
        @test mstr.table[2, 3] == idx

        @test mstr.table[1, 3] == 0
        @test mstr.table[1, 4] == 0
        StarAlgebras.cache!(mstr, [1], [3, 4])
        @test mstr.table[1, 3] == 3
        @test mstr.table[1, 4] == 4

        @test_throws StarAlgebras.ProductNotDefined mstr[k+1, k]
    end

    tmstr = StarAlgebras.CachedMTable{true}(b, table_size = (k, k))

    @test all(iszero, tmstr.table)
    @test tmstr[1, 2] == 2
    @test tmstr[2, 3] == b[star(b[2])*b[3]]
    @test tmstr[3, 2] == b[star(b[3])*b[2]]
end
