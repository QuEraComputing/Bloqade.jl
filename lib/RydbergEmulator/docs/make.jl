using Documenter
using RydbergEmulator
using DocThemeIndigo

indigo = DocThemeIndigo.install(RydbergEmulator)

format = if "latex" in ARGS
    Documenter.Writers.LaTeXWriter.LaTeX(platform="docker")
else
    Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://Happy-Diode.github.io/RydbergEmulator.jl",
        assets=String[indigo],
    )
end

makedocs(;
    modules = [RydbergEmulator],
    format,
    pages = [
        "Home" => "index.md",
        "Quick Start" => "quick-start.md",
        "References" => "ref.md",
    ],
    repo="https://github.com/Happy-Diode/RydbergEmulator.jl",
    sitename = "RydbergEmulator.jl",
)

deploydocs(;
    repo="github.com/Happy-Diode/RydbergEmulator",
)
