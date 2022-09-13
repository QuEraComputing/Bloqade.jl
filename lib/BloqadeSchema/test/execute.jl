using Test
using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema
using JSON
using Configurations


@testset "to_schema" begin
    T = 2
    atoms = [(π*i,0) for i in 1:10]

    # values = [1.0,nothing,Waveform(t->sin(2π*t/T)^2,T/2)]

    wf1 = Waveform(t->sin(π*t/T)^2,T)
    wf2 = Waveform(t->cos(π*t/T),T)
    wfs_mixed = [(i%2==0 ? 1 : 0.1)*wf1 + wf2 for i in 1:length(atoms)]
    wfs_same = [i*wf1+wf2 for i in 1:length(atoms)]
    scalar_values = [constant(;duration=T,value=1.0),wf1]
    all_values = [constant(;duration=T,value=1.0),wf1,wfs_mixed,wfs_same]
    params = get_device_capabilities_SI()

    check_atom_res(x) = all(BloqadeSchema.check_resolution.(params.lattice.geometry.positionResolution,x))
    check_clock_res(x) = all(BloqadeSchema.check_resolution.(params.rydberg.global_value.timeResolution,x))
    check_Δ_res(x) = all(BloqadeSchema.check_resolution.(params.rydberg.global_value.detuningResolution,x))
    check_δ_res(x) = all(BloqadeSchema.check_resolution.(params.rydberg.local_value.commonDetuningResolution,x))
    check_Δi_res(x) = all(BloqadeSchema.check_resolution.(params.rydberg.local_value.localDetuningResolution,x))
    check_ϕ_res(x) = all(BloqadeSchema.check_resolution.(params.rydberg.global_value.phaseResolution,x))
    check_Ω_res(x) = all(BloqadeSchema.check_resolution.(params.rydberg.global_value.rabiFrequencyResolution,x))

    for Ω in scalar_values, ϕ in scalar_values, Δ in all_values
        H = rydberg_h(atoms;Ω=Ω,Δ=Δ,ϕ=ϕ)
        j = BloqadeSchema.to_json(H,warn=true)
        t = Configurations.from_dict(BloqadeSchema.TaskSpecification, JSON.parse(j))

        sites = [site for (i,site) in enumerate(t.lattice.sites) if t.lattice.filling[i] == 1]
        
        rabi_freq_amp = t.effective_hamiltonian.rydberg.rabi_frequency_amplitude.global_value
        rabi_freq_phase = t.effective_hamiltonian.rydberg.rabi_frequency_phase.global_value
        detuning_global = t.effective_hamiltonian.rydberg.detuning.global_value
        detuning_local = t.effective_hamiltonian.rydberg.detuning.local_value

        @test all(check_atom_res(s) for s in sites)
        @test check_clock_res(rabi_freq_phase.times)
        @test check_clock_res(rabi_freq_amp.times)
        @test check_clock_res(detuning_global.times)
        @test check_ϕ_res(rabi_freq_phase.values)
        @test check_Ω_res(rabi_freq_amp.values)
        @test check_Δ_res(detuning_global.values)

        if !isnothing(detuning_local)
            @test check_clock_res(detuning_local.times)
            @test check_δ_res(detuning_local.values)
            @test check_Δi_res(detuning_local.lattice_site_coefficients)                  
        end


    end

end


# @testset "execute" begin
#     # Test if code will run without failing
#     Ω = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[0,1,1,0])
#     Δ = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[1,1,-1,-1])
#     ϕ = constant(;duration=3,value=0)
#     atoms = 5.0 * [i for i in 1:10]
    
#     H = rydberg_h(atoms,Ω=Ω,Δ=Δ,ϕ=ϕ)
#     task_string = to_json(H,n_shots=10,warn=true)
#     task_dict = BloqadeSchema.to_dict(H,n_shots=10,warn=true)
#     task = BloqadeSchema.to_schema(H,n_shots=10,warn=true)
    
#     r_string = execute(task_string)
#     r_dict = execute(task_dict)
#     r_task = execute(task)

# end
