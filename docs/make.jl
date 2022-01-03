using EaRyd
using EaRydCore
using EaRydODEEvolution
using Documenter
using DocThemeIndigo

indigo = DocThemeIndigo.install(EaRyd)
DocMeta.setdocmeta!(EaRyd, :DocTestSetup, :(using EaRyd); recursive=true)

makedocs(;
    modules=[EaRyd, EaRydCore, EaRydODEEvolution],
    authors="QuEra Computing Inc.",
    repo="https://github.com/Happy-Diode/EaRyd.jl/blob/{commit}{path}#{line}",
    sitename="EaRyd.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Happy-Diode.github.io/EaRyd.jl",
        assets=String[indigo],
    ),
    pages=[
        "Home" => "index.md",
        "Quick Start" => "quick-start.md",
        "The Julia Programming Language" => "julia.md",
        "Setting Up Company Registry" => "registry.md",
        "Tutorials" => [
            "Solving Maximum-Independent Set Using Rydberg QAOA" => "tutorials/mis.md"
        ],
        "References" => "ref.md",
    ],
)

deploydocs(;
    repo="github.com/Happy-Diode/EaRyd.jl",
)
