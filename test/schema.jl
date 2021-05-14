using Random
using Test
using Configurations
using RydbergEmulator

job = PulseJob(;
    natoms=10,
    radius=1.0,
    ff=0.8,
    pulses=[
        Pulse(
            duration=8e-4,
            rydberg=RydbergPulse(
                rabi=2π*4,
                detuning=0.2,
            )
        )
    ]
);

Random.seed!(1234)
atoms = square_lattice(10, 0.8)
h = RydInteract(atoms) + XTerm(10, 2π * 4) + NTerm(10, 0.2)
space = blockade_subspace(atoms, 1.0)
@test emulate(space, [8e-4], [h]) ≈ emulate(job) 

@test_throws AssertionError PulseJob(;
    lattice_constant_a = 0.1,
    lattice_constant_b = 10,
    pulses=[]
)

@test_throws AssertionError PulseJob(;
    lattice_constant_a = 10,
    lattice_constant_b = 2,
    pulses=[]
)

@test_throws AssertionError PulseJob(;
    lattice_angle = 240,
    lattice_constant_a = 10,
    lattice_constant_b = 20,
    pulses=[]
)

@test_throws AssertionError PulseJob(;
    natoms = 240,
    lattice_constant_a = 10,
    lattice_constant_b = 20,
    pulses=[]
)
