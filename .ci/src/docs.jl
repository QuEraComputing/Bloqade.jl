"""
documentation commands
"""
@cast module Doc

using Pkg
using Comonicon
using LiveServer
using ..BloqadeCI: root_dir, dev
using ..BloqadeCI.Example

function tutorial_pages()
    tutorials = Pair{String,String}[]
    Example.foreach_example() do path
        name = basename(path)
        mainjl = joinpath(path, "main.jl")
        target = joinpath("tutorials", name, "main.md")
        firstline = readline(mainjl)
        startswith(firstline, "# # ") || error("expecting example script start with # # <title>")
        return push!(tutorials, firstline[5:end] => target)
    end
    return tutorials
end

function pages(; light = false)
    PAGES = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "The Julia Programming Language" => "julia.md",
        "Manual" => [
            "Lattices" => "lattices.md",
            "Waveforms" => "waveform.md",
            "Hamiltonians" => "hamiltonians.md",
            "Registers and Observables" => "registers.md",
            "Emulation" => "emulation.md",
            "Working with Subspace" => "subspace.md",
            "Working with Units" => "units.md",
            "Maximum Independent Set" => "mis.md",
            "GPU Acceleration" => "cuda.md",
            "3-Level Support and Quantum Gates" => "3-level.md",
            "Interacting with Neutral Atom Hardware" => "schema.md",
            "Hardware Capabilities" => "capabilities.md"
        ],
    ]

    light || push!(PAGES, "Tutorials" => tutorial_pages())

    append!(PAGES, ["Contributing to Bloqade" => "contrib.md"])

    return PAGES
end

function render_all_examples()
    return Example.buildall(build_dir = root_dir("docs", "src", "tutorials"), target = "markdown", eval = true)
end

function doc_build_script(pages, repo)
    yao_pkgs = ["Yao", "YaoAPI", "YaoBlocks", "YaoArrayRegister", "Unitful"]
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
            canonical="https://QuEraComputing.github.io/Bloqade.jl",
            assets=String[indigo, "assets/favicon.ico"],
            sidebar_sitename=false,
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
    build_script = doc_build_script(pages(; light), "QuEraComputing/Bloqade.jl")
    return write(root_dir("docs", "make.jl"), build_script)
end

function setup_docs(light)
    dev("docs")
    light || dev_examples()
    generate_makejl(light)
    light || render_all_examples()
    return
end

@cast function serve(; host::String = "0.0.0.0", port::Int = 8000, light::Bool = false)
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
@cast function build(; light::Bool = false)
    setup_docs(light)
    docs_dir = root_dir("docs")
    docs_make_jl = root_dir("docs", "make.jl")
    julia_cmd = "using Pkg; Pkg.instantiate()"
    run(`$(Base.julia_exename()) --project=$docs_dir -e $julia_cmd`)
    return run(`$(Base.julia_exename()) --project=$docs_dir $docs_make_jl`)
end

end
