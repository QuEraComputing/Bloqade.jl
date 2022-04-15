function to_hamiltonian(task::TaskSpecification)
end

function to_json(h::AbstractBlock)
    # 1. check if the input block expression is a summation
    # of RydInteract, SumOfX, SumOfXPhase, SumOfN
    # 2. extract the atom positions, Ω, ϕ, Δ
    # 3. generate corresponding pulses in QuEraSchema
    # 4. call to_dict on QuEraSchema then using JSON3
    #    to convert the dict to json string
    # NOTE: let's ignore local pulse pattern for now
    # since in our representation it's all local pulse
end
