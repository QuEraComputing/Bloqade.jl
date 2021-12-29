# TODO: use GarishPrint after it gets smarter
function Base.show(io::IO, mime::MIME"text/plain", prob::DiscreteEvolution{P}) where P
    indent = get(io, :indent, 0)
    tab(indent) = " "^indent
    println(io, tab(indent), "DiscreteEvolution{", P, "}:")
    
    # state info
    print(io, tab(indent), "  reg: ")
    printstyled(io, typeof(prob.reg); color=:green)
    println(io)
    print(io, tab(indent), "  reg storage: ")
    printstyled(io, Base.format_bytes(sizeof(Yao.state(prob.reg))); color=:yellow)
    println(io)
    println(io)

    print(io, tab(indent), "  total duration: ")
    println(io, sum(prob.durations), " μs")

    println(io)
    println(io, tab(indent), "  hamiltonian: ")
    if length(prob.hs) < 6
        for h in prob.hs
            show(IOContext(io, :indent=>indent+4), mime, h)
            println(io)
        end
    else
        show(IOContext(io, :indent=>indent+4), mime, first(prob.hs))
        println(io)
        println(io)
        println(io, tab(indent), "       ⋮")
        println(io)
        show(IOContext(io, :indent=>indent+4), mime, last(prob.hs))
    end
    println(io)
    println(io)
    print(io, tab(indent), "  hamiltonian storage: ")
    printstyled(io, Base.format_bytes(storage_size(prob.cache)); color=:yellow)
    println(io)

    println(io)
    println(io, tab(indent), "  options: ")

    nfs = nfields(prob.options)
    for idx in 1:nfs
        name = fieldname(DiscreteOptions, idx)
        print(io, tab(indent), "    $name: ", repr(getfield(prob.options, idx)))
        if idx != nfs
            println(io)
        end
    end
end