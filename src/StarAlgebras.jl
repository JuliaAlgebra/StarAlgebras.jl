module StarAlgebras

using SparseArrays
import LinearAlgebra

import MutableArithmetics as MA

export StarAlgebra, AlgebraElement
export basis, coeffs, star

# AbstractCoefficients
## abstract definitions
include("coefficients.jl")
## concrete implementation
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

# Algebras and elts
include("types.jl")
include("algebra_elts.jl")
include("star.jl")

include("arithmetic.jl")
include("show.jl")

# augmented basis implementation
include("diracs_augmented.jl")

end
