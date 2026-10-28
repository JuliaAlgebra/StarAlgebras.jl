module TestParentPromotion

using Test, LinearAlgebra
import StarAlgebras as SA
import MutableArithmetics as MA

# A finite basis whose sorted elements can be extended without changing its type.
# Inserting smaller elements shifts keys, so promotion must reindex coefficients.
struct GrowingBasis <: SA.ExplicitBasis{Float64,Int}
    elements::Vector{Float64}
end

Base.length(b::GrowingBasis) = length(b.elements)
Base.getindex(b::GrowingBasis, i::Int) = b.elements[i]
Base.getindex(b::GrowingBasis, x::Float64) = searchsortedfirst(b.elements, x)
Base.:(==)(a::GrowingBasis, b::GrowingBasis) = a.elements == b.elements

function SA.zero_coeffs(::Type{T}, ::GrowingBasis) where {T}
    return SA.SparseCoefficients(Int[], T[])
end

function SA.promote_bases_with_maps(a::GrowingBasis, b::GrowingBasis)
    common = GrowingBasis(sort!(union(a.elements, b.elements)))
    map_a = Base.Fix1(getindex, common) ∘ Base.Fix1(getindex, a)
    map_b = Base.Fix1(getindex, common) ∘ Base.Fix1(getindex, b)
    return (common, map_a), (common, map_b)
end

function SA.promote_object(
    object::Float64,
    ::SA.DiracMStructure{Float64,Int,GrowingBasis},
    map,
)
    return object
end

function element(elements, values)
    alg = SA.StarAlgebra(1.0, GrowingBasis(elements))
    return SA.AlgebraElement(
        SA.SparseCoefficients(collect(eachindex(values)), values),
        alg,
    )
end

@testset "Dot promotes coefficient keys" begin
    a = element([2.0], ComplexF32[1+2im])
    b = element([1.0], ComplexF32[3+4im])
    c = element([1.0, 2.0], ComplexF32[0, 3+4im])
    alg_a, alg_b, alg_c = parent(a), parent(b), parent(c)
    # Key 1 denotes different basis elements in a and b.
    @test (@inferred dot(a, b)) === 0.0f0 + 0.0f0im
    @test (@inferred dot(a, c)) === 11.0f0 - 2.0f0im
    @test (@inferred dot(c, a)) === 11.0f0 + 2.0f0im
    @test (@inferred dot(a, a)) === 5.0f0 + 0.0f0im
    @test dot(zero(a), c) === 0.0f0 + 0.0f0im
    @test dot(a, zero(c)) === 0.0f0 + 0.0f0im
    @test parent(a) === alg_a
    @test parent(b) === alg_b
    @test parent(c) === alg_c
    @test collect(SA.nonzero_pairs(SA.coeffs(a))) == [1 => 1.0f0 + 2.0f0im]
    @test collect(SA.nonzero_pairs(SA.coeffs(b))) == [1 => 3.0f0 + 4.0f0im]
end

@testset "Mutable parent promotion" begin
    for T in (Int, BigInt)
        f = element([4.0, 8.0, 16.0], T[3, 1])
        a = element([2.0, 4.0, 8.0, 16.0], T[5])
        b = element([1.0, 2.0, 4.0, 8.0, 16.0], T[7, 2])
        original = MA.mutable_copy(f)
        alg = parent(f)
        t = SA.Term(alg, 1, T(3))
        @test parent(original) === alg
        @test parent(f + original) === alg

        for op in (+, -, *)
            output = MA.mutable_copy(f)
            @test (@inferred MA.operate!(op, output, a)) === output
            @test output == op(f, a)
            @test parent(output) != alg
            output = zero(f)
            @test (@inferred MA.operate_to!(output, op, a, b)) === output
            @test output == op(a, b)
        end
        for op in (+, -)
            output = zero(f)
            term = SA.Term(parent(a), 1, T(5))
            @test (@inferred MA.operate_to!(output, op, f, term)) === output
            @test output == op(f, a)
        end
        for op in (MA.add_mul, MA.sub_mul)
            output = zero(f)
            @test (@inferred MA.operate_to!(output, op, f, a, b, b)) === output
            @test output == op(f, a, b, b)
        end
        for op in (MA.add_mul, MA.sub_mul),
            left in (a, SA.Term(parent(a), 1, T(5)))

            output = MA.mutable_copy(f)
            @test (@inferred MA.operate!(op, output, left, b)) === output
            @test output == op(f, left, b)
            output = zero(f)
            @test (@inferred MA.operate_to!(output, op, f, left, b)) === output
            @test output == op(f, left, b)
            # The output can also be a factor while a later factor grows it.
            output = MA.mutable_copy(f)
            @test (@inferred MA.operate!(op, output, output, b)) === output
            @test output == op(f, f, b)
            output = MA.mutable_copy(f)
            alias = SA.AlgebraElement(SA.coeffs(output), parent(output))
            @test (@inferred MA.operate!(op, output, alias, b)) === output
            @test output == op(f, f, b)
        end
        output = MA.mutable_copy(f)
        @test (@inferred MA.add!!(output, SA.Term(parent(a), 1, T(5)))) ===
              output
        @test output == f + a
        @test f == original
        @test parent(f) === alg
        @test SA.basis(f).elements == [4.0, 8.0, 16.0]
        @test parent(t) === alg
        @test SA.basis_element(t) == 4.0
        @test SA.coeffs(a)[1] == 5
        @test SA.coeffs(b)[1] == 7
    end
end

# A basis with a type-level size exercises promotion that changes parent types.
struct PowerBasis{N} <: SA.ExplicitBasis{Float64,Int} end

Base.length(::PowerBasis{N}) where {N} = N
Base.getindex(::PowerBasis, i::Int) = 2.0^(i - 1)
Base.getindex(::PowerBasis, x::Float64) = Int(log2(x)) + 1
Base.:(==)(a::PowerBasis, b::PowerBasis) = typeof(a) === typeof(b)
function SA.zero_coeffs(::Type{T}, ::PowerBasis) where {T}
    return SA.SparseCoefficients(Int[], T[])
end

function SA.promote_bases_with_maps(
    a::PowerBasis{N},
    b::PowerBasis{M},
) where {N,M}
    if N == M
        return (a, nothing), (b, nothing)
    end
    common = PowerBasis{max(N, M)}()
    return (common, identity), (common, identity)
end

function SA.promote_object(
    object::Float64,
    ::SA.DiracMStructure{Float64,Int,<:PowerBasis},
    map,
)
    return object
end

const PowerAlgebra{N} = SA.StarAlgebra{
    Float64,
    Float64,
    SA.DiracMStructure{Float64,Int,PowerBasis{N},typeof(*)},
}
function Base.convert(::Type{PowerAlgebra{N}}, a::PowerAlgebra) where {N}
    return SA.StarAlgebra(SA.object(a), PowerBasis{N}())
end

@testset "Promotion between parent types" begin
    small = SA.StarAlgebra(1.0, PowerBasis{2}())
    large = SA.StarAlgebra(1.0, PowerBasis{4}())
    t, u = SA.Term(small, 2, 3), SA.Term(large, 3, 2.0)
    p, q = SA.algebra_element(t), SA.algebra_element(u)
    promoted_t = SA.Term(large, 2, 3.0)
    promoted_p = SA.algebra_element(promoted_t)
    @test (@inferred dot(p, q)) === 0.0
    @test (@inferred dot(p, promoted_p)) === 9.0
    @test (@inferred dot(t, u)) === 0.0
    @test (@inferred SA.promote_bases_with_maps(t, u)) ==
          ((promoted_t, identity), (u, identity))
    for (a, b, expected_a, expected_b) in (
        (t, u, promoted_t, u),
        (p, q, promoted_p, q),
        (p, u, promoted_p, q),
        (t, q, promoted_p, q),
    )
        promoted = @inferred promote(a, b)
        @test promoted == (expected_a, expected_b)
        @test typeof(promoted[1]) === typeof(expected_a)
        @test parent(promoted[1]) === large
        @test (@inferred convert(typeof(expected_a), a)) == expected_a
    end
    for op in (+, -, *)
        result = @inferred op(p, q)
        @test result == op(promoted_p, q)
        @test (@inferred MA.promote_operation(op, typeof(p), typeof(q))) ===
              typeof(result)
    end
    for x in (t, p), target in (large, SA.mstructure(large), SA.basis(large))
        promoted, common = @inferred SA.promote_bases(x, target)
        @test parent(promoted) === large
        @test promoted == (x isa SA.Term ? promoted_t : promoted_p)
        @test common == target
    end
    @test_throws InexactError convert(typeof(t), u)
    @test_throws InexactError convert(typeof(p), q)
    @test parent(t) === parent(p) === small
    @test SA.coeffs(p)[2] == 3
end

end # module
