# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module StarAlgebrasJuMPExt

import JuMP
import StarAlgebras

StarAlgebras.@scalar_arithmetic JuMP.AbstractJuMPScalar

end
