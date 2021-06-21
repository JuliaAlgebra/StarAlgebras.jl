# StarAlgebras

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kalmarek.github.io/StarAlgebras.jl/stable) -->
<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kalmarek.github.io/StarAlgebras.jl/dev) -->
[![Build Status](https://github.com/kalmarek/StarAlgebras.jl/workflows/CI/badge.svg)](https://github.com/kalmarek/StarAlgebras.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/kalmarek/StarAlgebras.jl/branch/main/graph/badge.svg?token=jpHVdYRx8G)](https://codecov.io/gh/kalmarek/StarAlgebras.jl)

----

The package implements `*`-algebras with basis. The prime example use is group/monoid algebras (or rings). An example usage can be as follows.

```julia
julia> using StarAlgebras

julia> using AbstractAlgebra

julia> G = SymmetricGroup(3)
Full symmetric group over 3 elements

julia> b = StarAlgebras.Basis{UInt8}(collect(G))
6-element Basis{Perm{Int64}, UInt8, Vector{Perm{Int64}}}:
 ()
 (1,2)
 (1,3,2)
 (2,3)
 (1,2,3)
 (1,3)

julia> RG = StarAlgebra(G, b)
*-Algebra of Full symmetric group over 3 elements

```

This creates the group algebra of the symmetric group. How do we compute inside the group algebra? There are a few ways to comstruct elements:
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
1·(1,3,2)

julia> f = AlgebraElement(rand(-3:3, order(Int, G)), RG) # an element given by vectors of coefficients in the basis
2·() -3·(1,2) +3·(1,3,2) -3·(2,3) -1·(1,2,3)

```
One may work with such element using the following functions:
```julia
julia> StarAlgebras.coeffs(f)
6-element Vector{Int64}:
  2
 -3
  3
 -3
 -1
  0

julia> StarAlgebras.star(p::Generic.Perm) = inv(p); star(f) # the star involution
2·() -3·(1,2) -1·(1,3,2) -3·(2,3) +3·(1,2,3)

julia> g = rand(G); StarAlgebras.coeffs(RG(g)) # note the type of coefficients
6-element SparseArrays.SparseVector{Int64, UInt8} with 1 stored entry:
  [1]  =  1

julia> x = one(RG) - 3RG(g); supp(x) # support of the funtion
1-element Vector{Perm{Int64}}:
 ()

julia> x(g) # value of x at g
-2

julia> x[g] += 3; x # modification of x in-place
1·()

julia> aug(f) # sum of coefficients
-2

julia> using LinearAlgebra; norm(f, 2) # Hilbert norm
5.656854249492381
```

Using this we can define e.g. a few projections in `RG` and check their orthogonality:
```julia
julia> using Test

julia> l = order(Int, G)
6

julia> P = sum(RG(g) for g in b) // l # projection to the subspace fixed by G
1//6·() +1//6·(1,2) +1//6·(1,3,2) +1//6·(2,3) +1//6·(1,2,3) +1//6·(1,3)

julia> @test P * P == P
Test Passed

julia> P3 = 2 * sum(RG(g) for g in b if sign(g) > 0) // l # projection to the subspace fixed by Alt(3) = C₃
1//3·() +1//3·(1,3,2) +1//3·(1,2,3)

julia> @test P3 * P3 == P3
Test Passed

julia> P2 = (RG(1) + RG(b[2])) // 2 # projection to the C₂-fixed subspace
1//2·() +1//2·(1,2)

julia> @test P2 * P2 == P2
Test Passed

julia> @test P2 * P3 == P3 * P2 == P # their intersection is precisely the same as the one for G
Test Passed

julia> P2m = (RG(1) - RG(b[2])) // 2 # orthogonal C₂-fixed subspace
1//2·() -1//2·(1,2)

julia> @test P2m * P2m == P2m
Test Passed

julia> @test iszero(P2m * P2) # indeed P2 and P2m are orthogonal
Test Passed

```


### More advanced use

`RG = StarAlgebra(G, b)` creates the algebra with `TrivialMStructure`, i.e. a multiplicative structure which computes product of basis elements every time it needs it. This of course may be wastefull, e.g. the computed products could be stored in a matrix for future use. There are two options here:
```julia
julia> mt = StarAlgebras.MTable(b, table_size=(length(b), length(b)))
6×6 StarAlgebras.MTable{UInt8, false, Matrix{UInt8}}:
 0x01  0x02  0x03  0x04  0x05  0x06
 0x02  0x01  0x04  0x03  0x06  0x05
 0x03  0x06  0x05  0x02  0x01  0x04
 0x04  0x05  0x06  0x01  0x02  0x03
 0x05  0x04  0x01  0x06  0x03  0x02
 0x06  0x03  0x02  0x05  0x04  0x01

```
creates an eagerly computed multiplication table on elements of `b`. Keyword `table_size` is used to specify the table size (above: it's the whole multiplication table). One can use the indexing syntax `mt[i,j]` to compute the **index** of the product of `i`-th and `j`-th elements of the basis. For example
```julia
julia> g = perm"(1,2,3)"; h = perm"(2,3)";

julia> i, j = b[g], b[h] # indices of g and h in basis b
(0x05, 0x04)

julia> k = mt[i,j] # the index of the product
0x06

julia> @test b[k] == g*h
Test Passed

```

The second option is
```julia
julia> cmt = StarAlgebras.CachedMTable(b, table_size=(length(b), length(b)));
```
This multiplication table is lazy, i.e. products will be computed and stored only when actually needed. Additionally, one may call
```julia
julia> using SparseArrays

julia> StarAlgebras.CachedMTable(b, spzeros(UInt8, length(b), length(b)));

```
to specify storage type of the matrix (by default it's a simple dense `Matrix`).
This may be advisable when a few products are computed repeatedly on a quite large basis.

```julia
julia> RGc = StarAlgebra(G, b, cmt)
*-Algebra of Full symmetric group over 3 elements

```
should be functinally equivalent to `RG` above, however it will cache computation of products lazily. A word of caution is needed here though. Even though `RGc` and `RG` are functionally equivalent, they are not **comparable** in the sense that e.g.
```julia
julia> @test one(RGc) != one(RG)
Test Passed

```

This is a conscious decision on our part, as comparing algebraic structures is easier said than done ;) To avoid solving this conundrum (are bases equal? are multiplicative structures equal? are these permuted by a compatible permutation? or maybe a linear transformation was applied to the basis, resulting in a different, but equivalent multiplicative structure?), elements could be mixed together **only if their parents are identically** (i.e. `===`) **equal**.

Finally, if the group is infinite (or just too large), but we need specific products, we may reduce the table_size to the required size (it doesn't have to be `length(b) × length(b)`). Note that in such case asking for a product outside of multiplication table will rise `ProductNotDefined` exception.

### Even more advanced use (for experts only)

This package originated as a tool to compute sum of (hermitian) squares in `*`-algebras. These consist not of standard `f*f` summands, but rather `star(f)*f`. You may think of semi-definite matrices: their Cholesky decomposition determines `P = Q'·Q`, where `Q'` denotes transpose. Algebra of matrices with transpose is an (the?) example of `*`-algebra.

To compute such sums of squares one may either sprinkle the code with `star`s, or define
```julia
julia> tcmt = StarAlgebras.CachedMTable{true}(b, table_size=(length(b), length(b)));

```
This multiplicative structure is **twisted** in the sense that `tcmt[i,j]` does not compute the product of `i`-th and `j`-th elements of the basis, but rather the `star` of `i`-th and `j`-th. An example
```julia
julia> k = tcmt[i,j]
0x02

julia> @test star(b[i])*b[j] == b[k]
Test Passed

```
you should only use it with extreme care, as this "product" is no longer associative! Observe:
```julia
julia> @test mt[mt[3, 5], 4] == mt[3, mt[5, 4]] # (b[3]*b[4])*b[5] == b[3]*(b[4]*b[5])
Test Passed

julia> @test tcmt[tcmt[3, 5], 4] == 0x06 # star(star(b[3])*b[5])*b[4] = star(b[5])*b[3]*b[4]
Test Passed

julia> @test tcmt[3, tcmt[5, 4]] == 0x04 # star(b[3])*star(b[5])*b[4]
Test Passed
```

However, writing sums of heritian squares is a breeze:
```julia
julia> tRG = StarAlgebra(G, b, tcmt)
*-Algebra of Full symmetric group over 3 elements

julia> x = tRG(perm"(1,2,3)")
1·(1,2,3)

julia> X = one(tRG) - x
1·() -1·(1,2,3)

julia> @test X^2 == X*X == 2one(tRG) - x - star(x)
Test Passed

```

#### WARNING!

Before using this mode you should consult the code (and potentially its author:) and understand very precisely what and where is happening!


-----
If you happen to use this package please cite either [1712.07167](https://arxiv.org/abs/1712.07167) or [1812.03456](https://arxiv.org/abs/1812.03456). This package superseeds [GroupRings.jl](https://github.com/kalmarek/GroupRings.jl) which was developed and used there. It served its purpose well. Let it rest peacefully.
