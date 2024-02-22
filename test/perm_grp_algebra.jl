using Pkg
Pkg.activate(@__DIR__)

using Revise
using StarAlgebras
import StarAlgebras as SA
using PermutationGroups
import PermutationGroups.AP as AP

SA.star(x::Number) = x'
SA.star(g::AP.AbstractPermutation) = inv(g)

G = PermGroup(perm"(1,2,3,4,5,6)", perm"(1,2)")
g, h = rand(G, 2)

db = SA.DiracBasis{UInt32}(G)
RG = SA.StarAlgebra(G, db)

ad = SA.AugmentedBasis(db)

ad[SA.AugmentedDirac(h)]

IG = SA.StarAlgebra(G, ad)

xcfs = SA.SparseCoefficients([one(G), g], [1, -1])
x = SA.AlgebraElement(xcfs, RG)

ycfs = SA.SparseCoefficients([one(G), inv(g)], [1, -1])
y = SA.AlgebraElement(ycfs, RG)

zcfs = SA.SparseCoefficients([one(G), h], [1, -1])
z = SA.AlgebraElement(zcfs, RG)

xycfs = SA.SparseCoefficients([one(G), g, inv(g)], [2, -1, -1])
xy = SA.AlgebraElement(xycfs, RG)

xzcfs = SA.SparseCoefficients([one(G), g, h, g * h], [1, -1, -1, 1])
xz = SA.AlgebraElement(xzcfs, RG)

@assert x != y
@assert x' == y
@assert SA.mstructure(basis(RG))(g, h) == SA.Dirac(g * h)
@assert x * y == xy
@assert x * z == xz

axcfs = SA.coeffs(coeffs(x), basis(RG), basis(IG))
aycfs = SA.coeffs(coeffs(y), basis(RG), basis(IG))
azcfs = SA.coeffs(coeffs(z), basis(RG), basis(IG))
ax = SA.AlgebraElement(axcfs, IG)
ay = SA.AlgebraElement(aycfs, IG)
az = SA.AlgebraElement(azcfs, IG)

@assert coeffs(ax * ay) == SA.coeffs(coeffs(x * y), basis(RG), basis(IG))
@assert coeffs(ax * az) == SA.coeffs(coeffs(x * z), basis(RG), basis(IG))

rcfs = SA.SparseCoefficients(rand(G, 100), rand(-2:2, 100))
r = SA.AlgebraElement(rcfs, RG)
scfs = SA.SparseCoefficients(rand(G, 100), rand(-2:2, 100))
s = SA.AlgebraElement(scfs, RG)

@time r * s

aug(r)
aug(s)
@assert aug(r * s) == aug(r) * aug(s)

m = PermutationGroups.order(UInt16, G)
fb = SA.FixedBasis(collect(G), SA.DiracMStructure(*), (m, m))
fRG = SA.StarAlgebra(G, fb)

fr = SA.AlgebraElement(coeffs(rcfs, basis(RG), basis(fRG)), fRG)
fs = SA.AlgebraElement(coeffs(scfs, basis(RG), basis(fRG)), fRG)

@assert aug(fr) == aug(r)
@assert aug(fs) == aug(s)

@assert aug(fr * fs) == aug(fr) * aug(fs)

r * s
fr * fs
frs = SA.AlgebraElement(coeffs(coeffs(r * s), basis(RG), basis(fRG)), fRG)
@assert fr * fs == frs

fr * fs
frs

let mt = SA.mstructure(basis(fRG)).table
    count(i -> isassigned(mt, i), eachindex(mt)), length(mt)
end
using BenchmarkTools
@btime $(star(r)) * $s
@btime $(star(fr)) * $fs

star(star(fr))

@btime star($fr)
@btime star($r)
@edit star(fr)

mulN(a, b, N) = [a * b for i in 1:N]

# using Profile

@profview mulN(fr, fs, 1000)
