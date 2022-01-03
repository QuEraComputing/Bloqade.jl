using Documenter
using EaRydCore
using DocThemeIndigo

indigo = DocThemeIndigo.install(EaRydCore)

format = if "latex" in ARGS
    Documenter.Writers.LaTeXWriter.LaTeX(platform="docker")
else
    Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://Happy-Diode.github.io/EaRydCore.jl",
        assets=String[indigo],
    )
end

makedocs(;
    modules = [EaRydCore],
    format,
    pages = [
        "Home" => "index.md",
        "Quick Start" => "quick-start.md",
        "References" => "ref.md",
    ],
    repo="https://github.com/Happy-Diode/EaRydCore.jl",
    sitename = "EaRydCore.jl",
)

deploydocs(;
    repo="github.com/Happy-Diode/EaRydCore",
)
