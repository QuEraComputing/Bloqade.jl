module BloqadeSchema

using Configurations

export TaskSpecification

include("types.jl")
include("serialize.jl")
include("deserialize.jl")

end
