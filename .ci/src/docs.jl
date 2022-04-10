"""
documentation commands
"""
@cast module Doc

using Pkg
using Comonicon
using LiveServer
using ..EaRydCI: root_dir, dev

function dev_examples()
    @info "setting up example Manifest.toml to use local packages"
    for each in readdir()
        project_dir = root_dir("examples", each)
        isdir(project_dir) || continue
        dev(project_dir)
    end
end

@cast function serve(;host::String="0.0.0.0", port::Int=8000)
    # setup environment
    dev("docs")
    dev_examples()
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
@cast function build()
    dev("docs")
    dev_examples()
    docs_dir = root_dir("docs")
    docs_make_jl = root_dir("docs", "make.jl")
    julia_cmd = "using Pkg; Pkg.instantiate()"
    
    for each in readdir(root_dir("examples"))
        example_dir = root_dir("examples", each)
        isdir(example_dir) || continue
        dev(example_dir)
    end
    run(`$(Base.julia_exename()) --project=$docs_dir -e $julia_cmd`)
    run(`$(Base.julia_exename()) --project=$docs_dir $docs_make_jl`)
end

end
