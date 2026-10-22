# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module StarAlgebrasMathOptInterfaceExt

import MathOptInterface as MOI
import StarAlgebras as SA

SA.@scalar_arithmetic MOI.AbstractScalarFunction

const ObjectWithoutIndex = Union{
    SA.Term{<:MOI.Utilities.ObjectWithoutIndex},
    SA.AlgebraElement{<:MOI.Utilities.ObjectWithoutIndex},
    SA.AbstractBasis,
}
# Empty tuples are already handled by MOI, independently of this extension.
const ObjectOrTupleWithoutIndex = Union{
    ObjectWithoutIndex,
    Tuple{ObjectWithoutIndex,Vararg{ObjectWithoutIndex}},
}
const ObjectOrTupleOrArrayWithoutIndex =
    Union{ObjectOrTupleWithoutIndex,AbstractArray{<:ObjectOrTupleWithoutIndex}}

function MOI.Utilities.map_indices(
    ::F,
    x::ObjectOrTupleOrArrayWithoutIndex,
) where {F<:Function}
    return x
end

function MOI.Utilities.substitute_variables(
    ::F,
    x::ObjectOrTupleOrArrayWithoutIndex,
) where {F<:Function}
    return x
end

end
