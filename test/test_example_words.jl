struct Word{T}
    alphabet::Vector{T}
    letters::Vector{Int}

    function Word(a::AbstractVector{T}, l::AbstractVector{<:Integer}) where {T}
        all(i -> 1 <= i <= length(a), l) ||
            throw(ArgumentError("Invalid word over alphabet $a: $w"))
        return new{T}(a, l)
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
    # star(:a) = :b
    # star(:b) = :a
    # star(:c) = :c

    star_d = Dict(1 => 2, 2 => 1)

    newletters = [get(star_d, l, l) for l in Iterators.reverse(w.letters)]
    return Word(w.alphabet, newletters)
end

function words(alphabet; radius)

    words = [Word(alphabet, Int[])] # word identity

    for r in 1:radius
        append!(
            words,
            [
                Word(alphabet, collect(w)) for
                w in Iterators.product(fill(1:length(alphabet), r)...)
            ],
        )
    end

    return words
end
