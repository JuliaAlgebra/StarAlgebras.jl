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

function Base.:(==)(w::Word, v::Word)
    return w.alphabet == v.alphabet && w.letters == v.letters
end
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

Base.isless(w::Word, v::Word) = w.letters < v.letters

struct FreeWords{T}
    alphabet::Vector{T}
end

Base.one(fw::FreeWords) = Word(fw.alphabet, Int[])

Base.eltype(::Type{FreeWords{T}}) where {T} = Word{T}
Base.IteratorSize(::Type{<:FreeWords}) = Base.IsInfinite()
function Base.iterate(aw::FreeWords)
    w = Word(aw.alphabet, Int[])
    stack = [w]
    return w, (stack, 1)
end

function Base.iterate(aw::FreeWords, state)
    stack, l = state
    if l > length(aw.alphabet)
        popfirst!(stack)
        l = 1
    end
    w = first(stack)
    nw = Word(aw.alphabet, [w.letters; l])
    push!(stack, nw)
    return nw, (stack, l + 1)
end

nwords(M::FreeWords, maxl::Integer) = nwords(M, 0, maxl)
function nwords(M::FreeWords, minl::Integer, maxl::Integer)
    maxl < minl && return zero(maxl)
    k = oftype(maxl, length(M.alphabet))
    return sum(k^i for i = minl:maxl)
end
