# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 10:37:33 2022
"""
# Python code to compare to
# C:\Users\jwurtz.QUERA\OneDrive - QuEra Computing\Documents\EaRyd\EaRyd.jl\examples\mis_adiabatic\main.jl

from numpy import *
from matplotlib.pyplot import *
import itertools

import quspin # https://weinbe58.github.io/QuSpin/
from pulser.devices import MockDevice




# Define locations of atoms as a ring
ground_state_spacing = 1    # How many nearest neighbors the unit disk touches.
                            # 1 is nearest neighbor, 2 is next nearest, etc.
                            # The ground state should be a Z(#+1) state.
R_interatomic = 9           # The unit disk radius, in um
L = 14                      # Number of atoms

# Compute the radius and positions of the atoms
radius = R_interatomic/(2*sin(pi/L * (ground_state_spacing)))
th = linspace(0,2*pi,L+1)[0:-1]
positions = array([radius * sin(th) , radius * cos(th)]).T

print('\nParameters of the atoms:')
print('Number of atoms:  {:0.0f}'.format(L))
print('Distance of atoms from origin: {:0.6f}um'.format(radius))
print('Unit disk radius: {:0.6f}um'.format(R_interatomic))

#%%
#positions = array([(0,0), (0,1), (0,2)])
N = positions.shape[0]

# The Rydberg energy at the unit disk radius sets the maximum energy scale of delta and omega
rydberg_blockade_energy = MockDevice.interaction_coeff / R_interatomic**6

Omega_max = 0.9*rydberg_blockade_energy

delta_0 = -2*rydberg_blockade_energy
delta_f = 0.9*rydberg_blockade_energy




Trise     = 2 # Rise time, in msec
Tramp     = 6 # Ramp time, in msec
Tfall     = 2 # Fall time, in msec
Tmax      = Trise + Tramp + Tfall

# Construct operators

"""
nn operator. n = (1-z)/2
"""
# 1. Construct the distances between different positions
dist = sqrt((positions[:,0:1] - positions[:,0:1].T)**2 + (positions[:,1::] - positions[:,1::].T)**2)
dist[range(N),range(N)] = inf
# 2. Compute the Rydberg weights
rydberg_weight =  MockDevice.interaction_coeff/ dist**6 / 4 # Units of microns and microseconds
# Because z is offset from n by 1, there is a constant Z term which is needed to be added
offset = -sum(rydberg_weight,0)
offset0 = 0.5*sum(offset) # The zero point offset

# 3. Add static Z terms together...
z_offset_field = [[offset[i],i] for i in range(N)]
ZZ_field       = [[rydberg_weight[i,j], i, j] for i,j in itertools.product(range(N),range(N)) if i>j]
static = [ ["zz",ZZ_field] , ["z",z_offset_field] ]

# The constant offset is simply sum(offset)/2

"""
X operator.
"""
#omega_of_T = lambda t: sin(pi*t/Tmax)
omega1_of_T = lambda t:Omega_max*(t/Trise)
omega2_of_T = lambda t:Omega_max
omega3_of_T = lambda t:Omega_max*(1 - (t/Tfall))
x_field = [[1/2*function(i),i] for i in range(N)]

"""
Z operator. n = (1-z)/2
"""
#delta_of_T = lambda t: cos(pi*t/Tmax)
delta1_of_T = lambda t: delta_0
delta2_of_T = lambda t: delta_0 * (1 - t/Tramp) + delta_f * (t/Tramp)
delta3_of_T = lambda t: delta_f
z_field = [[1/2,i] for i in range(N)]

# Add dynamic fields together...
#dynamic = [ ["x",x_field,omega_of_T,[]] , ["z",z_field,delta_of_T,[]]]
dynamic1 = [ ["x",x_field,omega1_of_T,[]] , ["z",z_field,delta1_of_T,[]] ]
dynamic2 = [ ["x",x_field,omega2_of_T,[]] , ["z",z_field,delta2_of_T,[]] ]
dynamic3 = [ ["x",x_field,omega3_of_T,[]] , ["z",z_field,delta3_of_T,[]] ]




# Define the Hamiltonian object...
bass = quspin.basis.spin_basis_1d(N)
ham1 = quspin.operators.hamiltonian(static,dynamic1,N=N,basis=bass)
ham2 = quspin.operators.hamiltonian(static,dynamic2,N=N,basis=bass)
ham3 = quspin.operators.hamiltonian(static,dynamic3,N=N,basis=bass)


# Define the wavefunction as the all-zero state.
psi0 = zeros(bass.Ns)
psi0[bass.state_to_int('0'*N)] = 1

# Sanity check: the zero state does in fact have an energy of zero.
ham = quspin.operators.hamiltonian(static,[],N=N,basis=bass).tocsc()

# Ground state has an energy of zero
assert isclose(dot(conj(psi0),ham.dot(psi0)) , offset0 ),'Ground state does not have energy zero'
# All ones state has an energy equal to the sum over Rydberg states
assert isclose ( ham[bass.state_to_int('1'*N),bass.state_to_int('1'*N)] - offset0 , sum(4*rydberg_weight)/2 ),'Most excited state is wrong!'


#%%

# Time-evolve the wavefunction
t0 = time.time()
psi1 = ham1.evolve(psi0,0,array([Trise]))
psi2 = ham2.evolve(psi1[:,0],0,array([Tramp]))
psi3 = ham3.evolve(psi2[:,0],0,array([Tfall]))
t1 = time.time()
print('Time to compute evolved state: {:0.3f}sec'.format(t1-t0))
#%%

print(sort(abs(psi3)**2,0)[-10::])
#energy = real(dot(conj(psi1[:,0]),ham(Tmax).dot(psi1[:,0])) - sum(offset)/2)
#print('Evolved excitation energy: {:0.4f}'.format(energy-emin))