using BloqadeSchema
using Test

function has_docstring(obj, sig::Type = Union{})
    Docs.getdoc(obj, sig) === nothing || return true
    binding = Docs.Binding(parentmodule(obj), nameof(obj))
    for mod in Docs.modules
        dict = Docs.meta(mod)
        haskey(dict, binding) || continue
        multidoc = dict[binding]
        for msig in multidoc.order
            sig <: msig && return true
        end
    end
    return false
end



@testset "Docstring Present" begin
    for symbol in names(BloqadeSchema)

        obj = Base.eval(@__MODULE__, :(BloqadeSchema.$symbol))

        typeof(obj) <: Union{Module} && continue
        
        @testset "$symbol" begin
            @test has_docstring(obj)
        end
    end
end