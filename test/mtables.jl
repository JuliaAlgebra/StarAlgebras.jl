struct Word{T}
    alphabet::Vector{T}
    letters::Vector{Int}

    function Word(a::AbstractVector{T}, l::AbstractVector{<:Integer}) where T
        all(i->1<=i<=length(a), l) || throw(ArgumentError("Invalid word over alphabet $a: $w"))
        return new{T}(a,l)
    end
end
Base.length(w::Word) = length(w.letters)
function Base.show(io::IO, w::Word)
    if isone(w)
        print(io, "(id)")
    else
        join(io, w.alphabet[w.letters], "·")
    end
end

Base.:(==)(w::Word, v::Word) = w.alphabet == v.alphabet && w.letters == v.letters
Base.hash(w::Word, h::UInt) = hash(w.alphabet, hash(w.letters, hash(Word, h)))

Base.one(w::Word) = Word(w.alphabet, Int[])
Base.isone(w::Word) = length(w) == 0

function Base.:*(w::Word, z::Word)
    @assert w.alphabet == z.alphabet
    return Word(w.alphabet, [w.letters; z.letters])
end

function StarAlgebras.star(w::Word)
    newletters = similar(w.letters)
    A = w.alphabet

    # star(:a) = :b
    # star(:b) = :a
    # star(:c) = :c

    for (i, l) in enumerate(Iterators.reverse(w.letters))
        newletters[i] = if A[l] === :a
            2
        elseif A[l] === :b
            1
        elseif A[l] === :c
            3
        end
    end
    return Word(w.alphabet, newletters)
end

elts = let A = [:a, :b, :c]
    w = Word(A, [1,2,3,2,1])

    words = [Word(A, [1]), Word(A, [2]), Word(A, [3])]
    for r in 2:4
        append!(
            words,
            [Word(A, collect(w)) for w in Iterators.product(fill(1:3, r)...)]
        )
    end
    pushfirst!(words, one(first(words)))
    words
end

@testset "TrivialMStructure" begin
    b = StarAlgebras.Basis{UInt8}(elts)
    k = findfirst(w->length(w)==3, b)-1

    mstr = StarAlgebras.TrivialMStructure{false}(b)

    @test mstr[1, 1] == 1
    @test mstr[1, 2] == 2
    @test mstr[2, 1] == 2
    @test mstr[3, 1] == 3

    idx = b[b[2]*b[3]] # == 8
    @test b[2]*b[3] == b[idx]
    @test mstr[2, 3] == idx

    idx = b[b[3]*b[2]]
    @test mstr[3, 2] == idx

    tmstr = StarAlgebras.TrivialMStructure{true}(b)

end

