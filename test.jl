# write your EaRyd example with Literate.jl here


using EaRyd
using EaRyd.EaRydODE.OrdinaryDiffEq


# Construct a list of atom positions. This is a 12-atom unit disk graph
# Positions are scaled so that the UDG radius is 9 microns.
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

atoms = EaRydLattices.AtomList(positions)
EaRydLattices.img_atoms(atoms)

# Define fixed parameters of the adiabatic ramp
omega_max = 2.3 * 2*pi   # Omega (X term) frequency, in units of radians/microsecond (eg, MHz)
U         = omega_max / 2.3
delta_0   = -3*U         # Initial Delta (Z term) frequency, in units of radians/microsecond (eg, MHz)
delta_1   =  1*U         #  Final  Delta (Z term) frequency, in units of radians/microsecond (eg, MHz)
# Adiabatic ramp times
Trise     = 2.
Tsweep    = 2.
Tfall     = 2.
Tmax = Trise+Tsweep+Tfall

# The ramp is a set of linear ramps. First, the Omega ramps from 0 to its maximum, creating a transverse-Ising-like model
#  Then, the Delta term ramps from negative to positive, which inverts the ground state from the all-zero state to a MIS state.
#  Finally, Omega ramps from its maximum back to zero, leaving the state in a Z eigenstate (eg logical basis) to read out an MIS
omega = piecewise_linear(clocks=[0, Trise, Trise+Tsweep, Trise+Tsweep+Tfall ], values=[0.0, omega_max , omega_max , 0])
delta = piecewise_linear(clocks=[0, Trise, Trise+Tsweep, Trise+Tsweep+Tfall ], values=[ delta_0, delta_0, delta_1, delta_1])


# Define the time-dependent Hamiltonian
ham = rydberg_h(positions; C = 5420158.53,Δ = delta, Ω = omega)


# Define the initial wavefunction as the all-zero state
psi0 = zero_state(12)
odesolve = ODEEvolution(psi0, Tmax, ham; dt=1e-4, algo=AB3())

@time emulate!(odesolve)

# Look at the last few elements to compare the state to Pulser
sort(broadcast(abs,odesolve.reg.state[:,1]).^2)

#=
# We can do the same thing in the blockade subspace, which substantially reduces the Hilbert space size and thus speeds up the numerics
space = blockade_subspace(positions, 9.)
psi0_blockade = zero_state(space)

@time begin
    odesolve_blockade = ODEEvolution(psi0_blockade,Tmax,ham,progress=true,dt=1e-2)
    #emulate!(odesolve_blockade)
    data = []
    for info in odesolve_blockade
        rydberg_density = [real(expect(put(12, i=>Op.n), info.reg)) for i in 1:12]
        push!(data,sum(rydberg_density))
    end
end
plot(data)
=#
