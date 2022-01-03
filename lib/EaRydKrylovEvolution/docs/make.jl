using Documenter
using EaRydKrylovEvolution
using DocThemeIndigo

indigo = DocThemeIndigo.install(EaRydKrylovEvolution)

format = if "latex" in ARGS
    Documenter.Writers.LaTeXWriter.LaTeX(platform="docker")
else
    Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://Happy-Diode.github.io/EaRydKrylovEvolution.jl",
        assets=String[indigo],
    )
end

makedocs(;
    modules = [EaRydKrylovEvolution],
    format,
    pages = [
        "Home" => "index.md",
        "Quick Start" => "quick-start.md",
        "References" => "ref.md",
    ],
    repo="https://github.com/Happy-Diode/EaRydKrylovEvolution.jl",
    sitename = "EaRydKrylovEvolution.jl",
)

deploydocs(;
    repo="github.com/Happy-Diode/EaRydKrylovEvolution",
)
