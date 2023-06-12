using Bloqade
using BloqadeKrylov


# create problem Hamiltonian 
r = zero_state(5)
atoms = [(i, ) for i in 1:5]
h = rydberg_h(atoms; Ω=sin)

# * using CFET integrator 
#    algo options (naming convention CFET<N>_<s>)
#      N: order of the magnus expansion. 
#      s: number of exponential-time propogator
#     -> CFET2_1(), CFET4_2(), CFET6_5(), CFET8_11()
prob = CFETEvolution(r, 0.0:1e-2:0.1, h, algo=CFET4_2())

# run the emulation
emulate!(prob)