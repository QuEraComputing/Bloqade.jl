
function Base.show(io::IO, ::MIME"text/plain", h::Hamiltonian)
    tab(n) = " "^(n + get(io, :indent, 0))
    println(io, tab(0), "Hamiltonian")
    println(io, tab(2), "number of dynamic terms: ", count(x -> x !== Base.one, h.fs))
    print(io, tab(2), "storage size: ")
    return print(io, Base.format_bytes(sum(sizeof, h.ts)))
end
