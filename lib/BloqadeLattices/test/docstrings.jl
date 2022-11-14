using BloqadeLattices
using Test

@testset "Docstring Present" begin
    for symbol in names(BloqadeLattices)

        symbol ∈ [:AbstractLattice, :AbstractRegion] && continue

        obj = Base.eval(@__MODULE__, :(BloqadeLattices.$symbol))

        typeof(obj) <: Union{Module} && continue
            
        docs_md = Base.eval(@__MODULE__, :(@doc BloqadeLattices.$symbol))
        doc_str = docs_md.meta[:results]
        @testset "$symbol" begin
            @test length(doc_str) > 0
        end
    end
end