# This file is a part of StarAlgebras.jl. License is MIT: https://github.com/JuliaAlgebra/StarAlgebras.jl/blob/main/LICENSE
# Copyright (c) 2026: Marek Kaluba, Benoît Legat

module StarAlgebrasJuMPExt

import JuMP
import StarAlgebras as SA

SA.@scalar_arithmetic JuMP.AbstractJuMPScalar

function JuMP.value(f::Function, t::SA.Term{<:JuMP.AbstractJuMPScalar})
    return SA.Term(parent(t), t.index, JuMP.value(f, SA.coefficient(t)))
end

function JuMP.value(
    f::Function,
    p::SA.AlgebraElement{<:JuMP.AbstractJuMPScalar},
)
    value = Base.Fix1(JuMP.value, f)
    output = similar(p, Base.promote_op(value, eltype(p)))
    return SA.map_coefficients_to!(output, value, p)
end

function JuMP.value(
    p::Union{
        SA.Term{<:JuMP.AbstractJuMPScalar},
        SA.AlgebraElement{<:JuMP.AbstractJuMPScalar},
    },
)
    return JuMP.value(JuMP.value, p)
end

function JuMP.function_string(
    mime::MIME"text/plain",
    p::Union{SA.Term,SA.AlgebraElement},
)
    return sprint(show, mime, p)
end

function JuMP.function_string(
    mime::MIME"text/latex",
    p::Union{SA.Term,SA.AlgebraElement},
)
    return SA.trim_LaTeX(mime, sprint(show, mime, p))
end

end
