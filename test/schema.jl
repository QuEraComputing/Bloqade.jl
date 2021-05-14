using Test
using Configurations
using RydbergEmulator

job = PulseJob(;
    lattice_constant_b = 10,
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
