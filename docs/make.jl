using Documenter, RydbergEmulator

makedocs(;
    modules=[RydbergEmulator],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/Happy-Diode/RydbergEmulator.jl",
    sitename="RydbergEmulator.jl",
)

deploydocs(;
    repo="github.com/Happy-Diode/RydbergEmulator",
)
