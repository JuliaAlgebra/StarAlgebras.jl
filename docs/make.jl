using StarAlgebras
using Documenter

DocMeta.setdocmeta!(StarAlgebras, :DocTestSetup, :(using StarAlgebras); recursive=true)

makedocs(;
    modules=[StarAlgebras],
    authors="Marek Kaluba <kalmar@mailbox.org>",
    repo="https://github.com/kalmar@amu.edu.pl/StarAlgebras.jl/blob/{commit}{path}#{line}",
    sitename="StarAlgebras.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kalmar@amu.edu.pl.github.io/StarAlgebras.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kalmar@amu.edu.pl/StarAlgebras.jl",
)
