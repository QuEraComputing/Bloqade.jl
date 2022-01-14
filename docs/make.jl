using Pkg
using EaRyd
using EaRydCore
using EaRydODE
using EaRydLattices
using EaRydWaveforms
using EaRydPlots
using Documenter
using DocThemeIndigo
using Literate

for each in readdir(pkgdir(EaRyd, "examples"))
    project_dir = pkgdir(EaRyd, "examples", each)
    isdir(project_dir) || continue
    input_file = pkgdir(EaRyd, "examples", each, "main.jl")
    output_dir = pkgdir(EaRyd, "docs", "src", "tutorials")
    Pkg.activate(project_dir)
    Pkg.instantiate()
    @info "executing" input_file
    Literate.markdown(input_file, output_dir; name=each, execute=true)
end

Pkg.activate(@__DIR__)
indigo = DocThemeIndigo.install(EaRyd)
DocMeta.setdocmeta!(EaRyd, :DocTestSetup, :(using EaRyd); recursive=true)

makedocs(;
    modules=[EaRyd, EaRydCore, EaRydODE, EaRydLattices, EaRydWaveforms, EaRydPlots],
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
            "Quantum Scar" => "tutorials/quantum-scar.md",
            "Adiabatic Evolution" => "tutorials/adiabatic.md",
            "Quantum Approximate Optimization Algorithm" => "tutorials/qaoa.md",
            # "Solving Maximum-Independent Set Using Rydberg QAOA" => "tutorials/mis.md",
        ],
        "References" => "ref.md",
    ],
)

deploydocs(;
    repo="github.com/Happy-Diode/EaRyd.jl",
)
