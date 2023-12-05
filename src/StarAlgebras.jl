module StarAlgebras

using SparseArrays
import LinearAlgebra

export StarAlgebra, AlgebraElement
export aug, basis, coeffs, star, supp

include("bases.jl")

include("mstructures.jl")
include("mtables.jl")

include("coefficients.jl")

include("types.jl")
include("algebra_elts.jl")
include("arithmetic.jl")
include("show.jl")

end
