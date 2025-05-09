# StarAlgebras

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kalmarek.github.io/StarAlgebras.jl/stable) -->
<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kalmarek.github.io/StarAlgebras.jl/dev) -->
[![CI](https://github.com/JuliaAlgebra/StarAlgebras.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaAlgebra/StarAlgebras.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/JuliaAlgebra/StarAlgebras.jl/branch/main/graph/badge.svg?token=jpHVdYRx8G)](https://codecov.io/gh/JuliaAlgebra/StarAlgebras.jl)

----

The package implements `*`-algebras with basis. The prime example use are group
or monoid algebras (or their finite dimensional subspaces).
An example usage can be as follows.

```julia
julia> using StarAlgebras; import StarAlgebras as SA

julia> using PermutationGroups

julia> G = PermGroup(perm"(1,2)", perm"(1,2,3)")
Permutation group on 2 generators generated by
 (1,2)
 (1,2,3)

julia> b = SA.DiracBasis{UInt32}(G);

julia> RG = StarAlgebra(G, b)
*-algebra of PermGroup( (1,2), (1,2,3) )

```

This creates the group algebra of the symmetric group. There are a few ways to
comstruct its elements:
```julia
julia> zero(RG)
0·()

julia> one(RG) # the canonical unit
1·()

julia> RG(1) # the same
1·()

julia> RG(-5.0) # coerce a scalar to the ring
-5.0·()

julia> RG(rand(G)) # the indicator function on a random element of G
1·()

julia> f = sum(rand(-3:3)*RG(rand(G)) for _ in 1:12) # an element given by vectors of coefficients in the basis
-1·() +4·(1,2) +5·(1,2,3) -2·(1,3,2) -2·(1,3)

julia> SA.star(g::PermutationGroups.AbstractPermutation) = inv(g); star(f) # the star involution
-1·() +4·(1,2) +5·(1,3,2) -2·(1,2,3) -2·(1,3)

julia> f' # the same
-1·() +4·(1,2) +5·(1,3,2) -2·(1,2,3) -2·(1,3)

julia> f' == f
false

julia> f - f' # the alternating part
7·(1,2,3) -7·(1,3,2)

julia> StarAlgebras.aug(p::PermutationGroups.Permutation) = 1

julia> StarAlgebras.aug(f) # sum of coefficients
4

julia> StarAlgebras.aug(f - f') # sum of coefficients
0

julia> using LinearAlgebra; norm(f, 2) # 2-norm
25.0

```
By default `f` here is stored as (sparse) coefficients against `basis(RG)`:
```julia
julia> coeffs(f) # same as SA.coeffs(f, basis(RG))
StarAlgebras.SparseCoefficients{...}(Permutation{...}[(), (1,2), (1,2,3), (1,3,2), (1,3)], [-1, 4, 5, -2, -2])
```

Working directly with elements of `G` one can also write

```julia
julia> g = Permutation(perm"(1,2,3)", G)
(1,2,3)

julia> x = one(RG) - 3RG(g); supp(x) # support of the funtion
2-element Vector{Permutation{…}}:
 ()
 (1,2,3)

julia> x(g) # value of x at g
-3

julia> x(inv(g))
0
```

## Example: Projections in group algebra

Using these functions let's define a few projections in `RG` and check their
orthogonality:
```julia
julia> using Test

julia> l = length(basis(RG))
6

julia> P = sum(RG(g) for g in b) // l # projection to the subspace fixed by all elements of G
1//6·() +1//6·(2,3) +1//6·(1,2) +1//6·(1,2,3) +1//6·(1,3,2) +1//6·(1,3)

julia> @test P * P == P
Test Passed

julia> P3 = 2 * sum(RG(g) for g in b if sign(g) > 0) // l # projection to the subspace fixed by Alt(3) = C₃
1//3·() +1//3·(1,2,3) +1//3·(1,3,2)

julia> @test P3 * P3 == P3
Test Passed

julia> h = PermutationGroups.gens(G,1)
(1,2)

julia> P2 = (RG(one(G)) + RG(h)) // 2 # projection to the C₂-fixed subspace
1//2·() +1//2·(1,2)

julia> @test P2 * P2 == P2
Test Passed

julia> @test P2 * P3 == P3 * P2 == P # their intersection is precisely the same as the one for G
Test Passed

julia> P2m = (RG(1) - RG(h)) // 2 # projection onto the subspace orthogonal to C₂-fixed subspace
1//2·() -1//2·(1,2)

julia> @test P2m * P2m == P2m
Test Passed

julia> @test iszero(P2m * P2) # indeed P2 and P2m are orthogonal
Test Passed
```

This package originated as a tool to compute sum of hermitian squares in
`*`-algebras. These consist not of standard `f*f` summands, but rather
`star(f)*f`. You may think of semi-definite matrices: their Cholesky
decomposition determines `P = Q'·Q`, where `Q'` denotes (conjugate) transpose.
Algebra of matrices with transpose is an (the?) example of a `*`-algebra.
To compute such sums of squares one may either sprinkle the code with `star`s,
or `'` (aka `Base.adjoint` postfix symbol):
```julia
julia> x = RG(Permutation(perm"(1,2,3)", G))
1·(1,2,3)

julia> X = one(RG) - x
1·() -1·(1,2,3)

julia> X'
1·() -1·(1,3,2)

julia> X'*X
2·() -1·(1,2,3) -1·(1,3,2)

julia> @test X'*X == star(X)*X == 2one(X) - x - star(x)
Test Passed

```

## Fixed basis and translation between bases

`b = SA.DiracBasis{UInt32}(G)` takes an iterator (in this case a finite
permutation group `G`) and creates a basis with Dirac multiplicative structure.
This means a basis of a linear space with the multiplicative structure where
multiplication of two basis elements results in a third one. This computation
does not depend on the finiteness of `G`, i.e. is fully lazy.

When the basis is known a priori one can create an efficient `SA.FixedBasis`
which caches (memoizes) the results of multiplications. For example:

```julia
julia> fb = let l = length(basis(RG))
    StarAlgebras.FixedBasis(b; n=l, mt=l) # mt can be also smaller than l
end;

julia> fRG = StarAlgebra(G, fb)
*-algebra of PermGroup( (1,2), (1,2,3) )

```

Since the length of the basis is known this time, algebra elements can be stored simply as (sparse) vectors:

```julia
julia> coeffs(one(fRG))
6-element SparseArrays.SparseVector{Int64, Int64} with 1 stored entry:
  [1]  =  1

julia> AlgebraElement(ones(Int, length(basis(fRG)))//6, fRG)
1//6·() +1//6·(2,3) +1//6·(1,2) +1//6·(1,2,3) +1//6·(1,3,2) +1//6·(1,3)

```

To translate coefficients between bases one may call

```julia
julia> coeffs(P2, fb)
6-element SparseArrays.SparseVector{Rational{Int64}, Int64} with 2 stored entries:
  [1]  =  1//2
  [3]  =  1//2

julia> coeffs(coeffs(P2), b, fb) # same as above, from source basis to target basis
6-element SparseArrays.SparseVector{Rational{Int64}, Int64} with 2 stored entries:
  [1]  =  1//2
  [3]  =  1//2

```

Translation in the opposite direction is also possible

```julia
julia> fP2 = AlgebraElement(StarAlgebras.coeffs(P2,fb), fRG)
1//2·() +1//2·(1,2)

julia> StarAlgebras.coeffs(fP2, b)
StarAlgebras.SparseCoefficients{…}(Permutation{…}[(), (1,2)], Rational{Int64}[1//2, 1//2])

julia> P2_ = AlgebraElement(StarAlgebras.coeffs(fP2, b), RG)
1//2·() +1//2·(1,2)

julia> @test P2 == P2_
Test Passed

```

-----
If you happen to use this package please cite either [1712.07167](https://arxiv.org/abs/1712.07167) or [1812.03456](https://arxiv.org/abs/1812.03456). This package superseeds [GroupRings.jl](https://github.com/kalmarek/GroupRings.jl) which was developed and used there. It served its purpose well. Let it rest peacefully.
