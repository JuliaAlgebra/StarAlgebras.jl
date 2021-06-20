module StarAlgebras

using SparseArrays
import LinearAlgebra

include("bases.jl")

include("mstructures.jl")
include("mtables.jl")

include("types.jl")
include("algebra_elts.jl")
include("arithmetic.jl")
include("show.jl")

end
