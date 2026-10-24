# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module StarAlgebrasGroupsCoreExt

import GroupsCore
import StarAlgebras as SA

SA.@object_arithmetic GroupsCore.MonoidElement

end
