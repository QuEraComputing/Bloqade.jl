module BloqadePython

using Bloqade
using PythonCall

include("interface.jl")

function __init__()
    init_jlwrap_block()
end

end
