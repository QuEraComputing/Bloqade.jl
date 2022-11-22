using Test
using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema
using JSON
using Configurations
using Unitful: μs, s, MHz, rad 
using Logging


@testset "ViolationException" begin
    
    Ω = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[0,1,1,0])
    Δ = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[1,1,-1,-1])
    ϕ = piecewise_constant(;clocks=Float64[0,1,2,4],values=Float64[0,1,-1])
    atoms = 5.0 * [i for i in 1:10]
    h = rydberg_h(atoms, Ω=Ω, Δ=Δ, ϕ=ϕ)
    @test_throws BloqadeSchema.ValidationException to_json(h)

end

@testset "convert_units" begin
    T = 1 # seconds
    T_list = [i for i in 1:10]
    pair = (1,2)
    pair_list = [(i,j) for i in 1:5 for j in 1:5]
    @test 1e6 ≈ BloqadeSchema.convert_units(T,s,μs)
    @test (1e6 .* T_list) ≈ BloqadeSchema.convert_units(T_list,s,μs)
    @test all((1e6 .* pair ).≈ BloqadeSchema.convert_units.(pair,s,μs))

end

@testset "to_schema"  begin
    
    T = 0.5
    atoms = [(5*i,0) for i in 1:10]

    # values = [1.0,nothing,Waveform(t->sin(2π*t/T)^2,T/2)]

    wf1 = Waveform(t->sin(π*t/T)^2,T)
    wf2 = Waveform(t->cos(π*t/T),T)
    wfs_mixed = [(i%2==0 ? 1 : 0.1)*wf1 + wf2 for i in 1:length(atoms)]
    wfs_same = [i*wf1+wf2 for i in 1:length(atoms)]
    scalar_values = [constant(;duration=T,value=1.0),wf1,wf2]
    # all_values = [constant(;duration=T,value=1.0),wf1,wfs_mixed,wfs_same]
    params = get_device_capabilities_SI()

    check_atom_res(x) = !any(BloqadeSchema.check_resolution.(params.lattice.geometry.position_resolution,x))
    check_clock_res(x) = !any(BloqadeSchema.check_resolution.(params.rydberg.global_value.time_resolution,x))
    check_Δ_res(x) = !any(BloqadeSchema.check_resolution.(params.rydberg.global_value.detuning_resolution,x))
    # check_δ_res(x) = !any(BloqadeSchema.check_resolution.(params.rydberg.local_value.common_detuning_resolution,x))
    # check_Δi_res(x) = !any(BloqadeSchema.check_resolution.(params.rydberg.local_value.local_detuning_resolution,x))
    check_ϕ_res(x) = !any(BloqadeSchema.check_resolution.(params.rydberg.global_value.phase_resolution,x))
    check_Ω_res(x) = !any(BloqadeSchema.check_resolution.(params.rydberg.global_value.rabi_frequency_resolution,x))

    for Ω in scalar_values, ϕ in scalar_values, Δ in scalar_values
        h = rydberg_h(atoms;Ω=Ω,Δ=Δ,ϕ=ϕ)
        h_hw,info = hardware_transform(h)
        j = BloqadeSchema.to_json(h_hw)

        t = Configurations.from_dict(BloqadeSchema.TaskSpecification, JSON.parse(j))

        sites = [site for (i,site) in enumerate(t.lattice.sites) if t.lattice.filling[i] == 1]
        
        rabi_freq_amp = t.effective_hamiltonian.rydberg.rabi_frequency_amplitude.global_value
        rabi_freq_phase = t.effective_hamiltonian.rydberg.rabi_frequency_phase.global_value
        detuning_global = t.effective_hamiltonian.rydberg.detuning.global_value
        detuning_local = t.effective_hamiltonian.rydberg.detuning.local_value

        pos_res = params.lattice.geometry.position_resolution
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


@testset "execute" begin
    # Test if code will run without failing
    Ω = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[0,1,1,0])
    Δ = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[1,1,-1,-1])
    ϕ = constant(;duration=3,value=0)
    atoms = 5.0 * [i for i in 1:10]
    
    h = rydberg_h(atoms,Ω=Ω,Δ=Δ,ϕ=ϕ)
    h,info = hardware_transform(h)
    task_string = to_json(h;n_shots=10)
    task_dict = BloqadeSchema.to_dict(h;n_shots=10)
    task = BloqadeSchema.to_schema(h;n_shots=10)
    
    @test execute(task_string) isa String
    @test execute(task_dict) isa String
    @test execute(task) isa String

end


@testset "conversion to and from schema" begin
        # Test if code will run without failing
        Ω = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[0,1,1,0])
        Δ = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[1,1,-1,-1])
        ϕ = constant(;duration=3,value=0)
        atoms = [(5.0*i,0.0) for i in 1:4]
        
        h = rydberg_h(atoms,Ω=Ω,Δ=Δ,ϕ=ϕ)
        h,info = hardware_transform(h)
        @test h == BloqadeSchema.from_schema(BloqadeSchema.to_schema(h))
        @test h == BloqadeSchema.from_dict(BloqadeSchema.to_dict(h))
        @test h == BloqadeSchema.from_json(BloqadeSchema.to_json(h))
end


@testset "to_json_no_validation" begin
    # Test if code will run without failing
    Ω = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[0,1,1,0])
    Δ = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[1,1,-1,-1])
    ϕ = piecewise_constant(;clocks=Float64[0,1,2,3],values=Float64[1,1,-1])
    atoms = [(5.0*i,0.0) for i in 1:4]
    
    h = rydberg_h(atoms,Ω=Ω,Δ=Δ,ϕ=ϕ)
    h,info = hardware_transform(h)    
     atoms,ϕ,Ω,Δ = get_rydberg_params(h)

    @test to_json_no_validation(atoms,ϕ=ϕ,Ω=Ω,Δ=Δ) == to_json(h)
    @test to_schema_no_validation(atoms,ϕ=ϕ,Ω=Ω,Δ=Δ) == to_schema(h)

    # test with partial filling
    coords = [(5.0*i,0.0) for i in 1:4]
    filling = [0,1,1,0]

    lattice = BloqadeSchema.Lattice(sites=atoms,filling=filling)

    @test to_json_no_validation(lattice,ϕ=ϕ,Ω=Ω,Δ=Δ) == "{\"nshots\":1,\"lattice\":{\"sites\":[[5.0e-6,0.0],[1.0e-5,0.0],[1.5e-5,0.0],[2.0e-5,0.0]],\"filling\":[0,1,1,0]},\"effective_hamiltonian\":{\"rydberg\":{\"rabi_frequency_amplitude\":{\"global\":{\"times\":[0.0,1.0e-6,2.0e-6,3.0e-6],\"values\":[0.0,1.0e6,1.0e6,0.0]}},\"rabi_frequency_phase\":{\"global\":{\"times\":[0.0,1.0e-6,2.0e-6,3.0e-6],\"values\":[0.0,1.0,-1.0,-1.0]}},\"detuning\":{\"global\":{\"times\":[0.0,1.0e-6,2.0e-6,3.0e-6],\"values\":[1.0e6,1.0e6,-1.0e6,-1.0e6]}}}}}"
end   

@testset "nshots validation" begin
    Ω = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[0,1,1,0])
    Δ = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[1,1,-1,-1])
    ϕ = constant(;duration=3,value=0)
    atoms = [(5.0*i,0.0) for i in 1:4]
    
    h = rydberg_h(atoms,Ω=Ω,Δ=Δ,ϕ=ϕ)
    h,info = hardware_transform(h)

    dc = get_device_capabilities()

    @test_throws BloqadeSchema.ValidationException to_json(h;n_shots=dc.task.number_shots_min-1)
    @test_throws BloqadeSchema.ValidationException to_json(h;n_shots=dc.task.number_shots_max+1)

end