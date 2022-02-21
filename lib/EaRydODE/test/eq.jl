using EaRydODE
using LinearAlgebra
using BenchmarkTools
using SparseArrays
using EaRydWaveforms
using EaRydCore
using EaRydODE: EquationCache

positions =
      [(  9.67604084,  -2.30465001),
       (  4.48561432,  -5.53855587),
       (  6.21928636,   1.67251608),
       ( -9.9279178 ,   0.67689399),
       (  5.73589863,   7.68914674),
       ( -6.10562758,  -4.61723791),
       ( -0.72276152, -12.47076581),
       (  0.4587603 ,   9.43256501),
       ( -4.89777907,   8.13329584),
       ( -0.99712562,  -5.45801548),
       ( -4.73810949,   1.78019481),
       (  0.81372062,   1.00461261)]

# Define fixed parameters of the adiabatic ramp
omega_max = 2.3 * 2*pi   # Omega (X term) frequency, in units of radians/microsecond (eg, MHz)
U         = omega_max / 2.3
delta_0   = -3*U         # Initial Delta (Z term) frequency, in units of radians/microsecond (eg, MHz)
delta_1   =  1*U         #  Final  Delta (Z term) frequency, in units of radians/microsecond (eg, MHz)
# Adiabatic ramp times
Trise     = 2
Tsweep    = 2
Tfall     = 2
Tmax = Trise+Tsweep+Tfall

# The ramp is a set of linear ramps. First, the Omega ramps from 0 to its maximum, creating a transverse-Ising-like model
#  Then, the Delta term ramps from negative to positive, which inverts the ground state from the all-zero state to a MIS state.
#  Finally, Omega ramps from its maximum back to zero, leaving the state in a Z eigenstate (eg logical basis) to read out an MIS
omega = piecewise_linear(clocks=[0, Trise, Trise+Tsweep, Trise+Tsweep+Tfall ], values=[0.0, omega_max , omega_max , 0])
delta = piecewise_linear(clocks=[0, Trise, Trise+Tsweep, Trise+Tsweep+Tfall ], values=[ delta_0, delta_0, delta_1, delta_1])


# Define the time-dependent Hamiltonian
ham = rydberg_h(positions; C = 5420158.53,Δ = delta, Ω = omega)

eq = SchrodingerEquation(ham, FullSpace(), EquationCache(SparseMatrixCSC(ham(0.1))))

dstate = zeros(ComplexF64, 1 << length(positions))
state = rand(ComplexF64, 1 << length(positions))
t = 0.1
eq.cache.hamiltonian
@benchmark $eq($dstate, $state, nothing, t)



tc = split_term(ComplexF64, ham, FullSpace())
function equation(tc, dstate, state, t)
    for (f, h) in zip(tc.fs, tc.hs)
        mul!(dstate, h, state, f(t), one(t))
    end
    lmul!(-im, dstate)
    return dstate
end

@benchmark equation($tc, $dstate, $state, 0.1)

dstate1 = zeros(ComplexF64, 1 << length(positions))
dstate2 = zeros(ComplexF64, 1 << length(positions))
eq(dstate1, state, nothing, t)
equation(tc, dstate2, state, 0.1)
dstate1 ≈ dstate2

put(5, 1 => X)
