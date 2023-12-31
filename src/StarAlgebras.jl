module StarAlgebras

using SparseArrays
import LinearAlgebra

import MutableArithmetics as MA

export StarAlgebra, AlgebraElement
export aug, basis, coeffs, star, supp

# AbstractCoefficients
## abstract definitions
include("coefficients.jl")
## concrete implementations
include("diracs.jl")
include("diracs_augmented.jl")
include("sparse_coeffs.jl")

# MultiplicativeStructures
include("mstructures.jl")
include("mtables.jl")

# AbstractBases
## abstract definitions
include("bases.jl")
# concrete implementations
include("bases_dirac.jl")
include("bases_fixed.jl")

# star depends only on basis and coefficients
include("star.jl")

# Algebras and elts
include("types.jl")
include("algebra_elts.jl")

include("arithmetic.jl")
include("show.jl")

end
