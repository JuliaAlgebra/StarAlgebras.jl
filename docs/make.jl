using StarAlgebras
using Documenter

DocMeta.setdocmeta!(StarAlgebras, :DocTestSetup, :(using StarAlgebras); recursive=true)

makedocs(;
    modules=[StarAlgebras],
    authors="Marek Kaluba and Benoît Legat",
    sitename="StarAlgebras.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaAlgebra/StarAlgebras.jl.git", push_preview=true, devbranch="main"
)
