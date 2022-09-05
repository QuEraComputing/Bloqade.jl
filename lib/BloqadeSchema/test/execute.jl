using Test
using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema
using JSON


@testset "to_schema" begin
    T = 2
    atoms = [(π*i,0) for i in 1:10]

    # values = [1.0,nothing,Waveform(t->sin(2π*t/T)^2,T/2)]
    values = [1.0,constant(;duration=T,value=1.0),Waveform(t->sin(2π*t/T)^2,T)]
    params = BloqadeSchema.SchemaConversionParams()

    check_atom_res = x->all(BloqadeSchema.check_resolution.(params.atom_position_resolution,x))
    check_clock_res = x->all(BloqadeSchema.check_resolution.(params.rabi_time_resolution,x))
    check_Δ_res = x->all(BloqadeSchema.check_resolution.(params.rabi_detuning_resolution,x))
    check_Δi_res = x->all(BloqadeSchema.check_resolution.(params.rabi_detuning_local_resolution,x))
    check_ϕ_res = x->all(BloqadeSchema.check_resolution.(params.rabi_frequency_phase_resolution,x))
    check_Ω_res = x->all(BloqadeSchema.check_resolution.(params.rabi_frequency_amplitude_resolution,x))

    for Ω in values
        for ϕ in values
            for Δ in values
                if any((f isa BloqadeSchema.DynamicParam) for f in [ϕ,Ω,Δ])
                    H = rydberg_h(atoms;Ω=Ω,Δ=Δ,ϕ=ϕ)
                    j = BloqadeSchema.to_json(H,waveform_tolerance=1e-2,warn=true,discretize=true)
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
                        @test check_Δ_res(detuning_local.values)
                        @test check_Δi_res(detuning_local.lattice_site_coefficients)                  
                    end


                end
            end
        end
    end

end

