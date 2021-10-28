using EaRyd
using Documenter

DocMeta.setdocmeta!(EaRyd, :DocTestSetup, :(using EaRyd); recursive=true)

makedocs(;
    modules=[EaRyd],
    authors="Roger-Luo <rogerluo.rl18@gmail.com> and contributors",
    repo="https://github.com/Happy-Diode/EaRyd.jl/blob/{commit}{path}#{line}",
    sitename="EaRyd.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Happy-Diode.github.io/EaRyd.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Happy-Diode/EaRyd.jl",
)
