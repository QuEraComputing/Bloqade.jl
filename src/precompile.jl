if ccall(:jl_generating_output, Cint, ()) == 1   # if we're precompiling the package
    let
        nsites = 4
        atoms = generate_sites(ChainLattice(), nsites, scale = 5.72)
        total_time = 3.0
        Ω_max = 2π * 4
        Ω = piecewise_linear(clocks = [0.0, 0.1, 2.1, 2.2, total_time], values = [0.0, Ω_max, Ω_max, 0, 0])
        # The detuning sequence can also be created in a similar way.
        U1 = -2π * 10
        U2 = 2π * 10
        Δ = piecewise_linear(clocks = [0.0, 0.6, 2.1, total_time], values = [U1, U1, U2, U2])
        show(devnull, MIME"text/plain"(), Δ)
        h = rydberg_h(atoms; Δ, Ω)
        show(devnull, MIME"text/plain"(), h)
        reg = zero_state(nsites)
        prob = SchrodingerProblem(reg, total_time, h)
        show(devnull, MIME"text/plain"(), prob)
        emulate!(prob)
    end
end
