using Test
using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema
using JSON


@testset "to_schema" begin
    T = 4
    atoms = [(i,0) for i in 1:10]

    values = [1.0,nothing,Waveform(t->sin(2π*t/T)^2,T/2)]
    for Ω in values
        for ϕ in values
            for Δ in values
                if any((f isa BloqadeSchema.DynamicParam) for f in [ϕ,Ω,Δ])
                    H = rydberg_h(atoms;Ω=Ω,Δ=Δ,ϕ=ϕ)
                    h = BloqadeSchema.to_json(H,waveform_tolerance=1e-2,warn=true)
                end
            end
        end
    end

end


@testset "execute" begin
    Ω = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[0,1,1,0])
    Δ = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[1,1,-1,-1])
    
    atoms = 5.0 * [i for i in 1:10]
    
    H = rydberg_h(atoms,Ω=Ω,Δ=Δ)
    task_string = to_json(H,n_shots=10)
    task_dict = BloqadeSchema.to_dict(H,n_shots=10)
    task = BloqadeSchema.to_schema(H,n_shots=10)
    
    r_string = execute(task_string)
    r_dict = execute(task_dict)
    r_task = execute(task)

end
