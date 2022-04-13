"""
documentation commands
"""
@cast module Doc

using Pkg
using Comonicon
using LiveServer
using ..BloqadeCI: root_dir, dev
using ..BloqadeCI.Example

const BEFORE_TUTORIAL = [
    "Home" => "index.md",
    "Installation" => "install.md",
    "The Julia Programming Language" => "julia.md",
    "Manual" => [
        "Waveforms" => "waveform.md",
        "Lattices" => "lattices.md",
        "Hamiltonians" => "hamiltonians.md",
        "Registers" => "registers.md",
        "Emulation" => "emulation.md",
        "Observables" => "observables.md",
        "Maximum Independent Set" => "mis.md",
        "CUDA Acceleration" => "cuda.md",    
    ],
]

const AFTER_TUTORIAL = [
    "Advanced Topics" => [
        "Rydberg Blockade" => "topics/blockade.md",
        "Bravais Lattice" => "topics/bravais.md",
        "Automatic Differentiation" => "topics/ad.md",
    ],
    "Contributing Bloqade" => "contrib.md",
]

const LIGHT_PAGES = [
    BEFORE_TUTORIAL...,
    AFTER_TUTORIAL...,
]

const PAGES=[
    BEFORE_TUTORIAL...,
    "Tutorials" => [
        "Quantum Scar" => "tutorials/quantum-scar/main.md",
        "Adiabatic Evolution" => "tutorials/adiabatic/main.md",
        "Quantum Approximate Optimization Algorithm" => "tutorials/qaoa/main.md",
        # "Solving Maximum-Independent Set Using Rydberg QAOA" => "tutorials/mis.md",
    ],
    AFTER_TUTORIAL...
]

function render_all_examples()
    Example.buildall(
        build_dir=root_dir("docs", "src", "tutorials"),
        target="markdown",
        eval=true,
    )
end

function doc_build_script(pages, repo)
    yao_pkgs = [
        "Yao", "YaoAPI", "YaoBase",
        "YaoBlocks", "YaoArrayRegister",
    ]
    using_stmts = ["Documenter", "DocThemeIndigo"]
    non_cuda_pkgs = filter(!isequal("BloqadeCUDA"), readdir(root_dir("lib")))
    push!(non_cuda_pkgs, "Bloqade")
    append!(non_cuda_pkgs, yao_pkgs)
    append!(using_stmts, non_cuda_pkgs)
    
    return """
    $(join(map(x->"using "*x, using_stmts), "\n"))

    indigo = DocThemeIndigo.install(Bloqade)
    DocMeta.setdocmeta!(Bloqade, :DocTestSetup, :(using Bloqade); recursive=true)

    makedocs(;
        root=\"$(root_dir("docs"))\",
        modules=[$(join(non_cuda_pkgs, ", "))],
        authors="QuEra Computing Inc.",
        repo="https://github.com/$repo/blob/{commit}{path}#{line}",
        sitename="Bloqade.jl",
        doctest=false,
        format=Documenter.HTML(;
            prettyurls=get(ENV, "CI", "false") == "true",
            canonical="https://Happy-Diode.github.io/Bloqade.jl",
            assets=String[indigo],
        ),
        pages=$pages,
    )
    deploydocs(;
        repo="github.com/$repo",
    )
    """
end

function dev_examples()
    @info "setting up example Manifest.toml to use local packages"
    for each in readdir()
        project_dir = root_dir("examples", each)
        isdir(project_dir) || continue
        dev(project_dir)
    end
end

function generate_makejl(light)
    build_script = if light
        doc_build_script(LIGHT_PAGES, "Happy-Diode/Bloqade.jl")
    else
        doc_build_script(PAGES, "Happy-Diode/Bloqade.jl")
    end

    write(root_dir("docs", "make.jl"), build_script)
end

function setup_docs(light)
    dev("docs")
    light || dev_examples()
    generate_makejl(light)
    light || render_all_examples()
    return
end

@cast function serve(;host::String="0.0.0.0", port::Int=8000, light::Bool=false)
    # setup environment
    setup_docs(light)

    docs_dir = root_dir("docs")
    julia_cmd = "using Pkg; Pkg.instantiate()"
    run(`$(Base.julia_exename()) --project=$docs_dir -e $julia_cmd`)

    docs_dir = root_dir(".ci")
    serve_cmd = """
    using LiveServer;
    LiveServer.servedocs(;
        doc_env=true,
        skip_dirs=[
            joinpath("docs", "src", "assets"),
            joinpath("docs", "src", "tutorials"),
        ],
        literate="examples",
        host=\"$host\",
        port=$port,
    )
    """
    try
        run(`$(Base.julia_exename()) --project=$(root_dir(".ci")) -e $serve_cmd`)
    catch e
        if e isa InterruptException
            return
        else
            rethrow(e)
        end
    end
    return
end

"""
build the docs
"""
@cast function build(;light::Bool=false)
    setup_docs(light)
    docs_dir = root_dir("docs")
    docs_make_jl = root_dir("docs", "make.jl")
    julia_cmd = "using Pkg; Pkg.instantiate()"
    run(`$(Base.julia_exename()) --project=$docs_dir -e $julia_cmd`)
    run(`$(Base.julia_exename()) --project=$docs_dir $docs_make_jl`)
end

end
