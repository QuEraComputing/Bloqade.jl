# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 09:06:07 2022

Copied and adapted from Pasqual pulser example
https://pulser.readthedocs.io/en/stable/tutorials/simulating.html

@author: jwurtz
"""
import qutip
qutip.settings.num_cpus = 1

from pulser import Pulse, Sequence, Register
from pulser.simulation import Simulation
from pulser.waveforms import BlackmanWaveform, RampWaveform
from pulser.devices import MockDevice

from numpy import *
import time



def pulser_compare(Nqubits , Ttotal ):
    """
    Simulates an adiabatic anneal on an Nqubit ring of Rydberg atoms,
    with each qubit seperated by 9 microns
    Nqubits - number of qubits
    Ttotal  - Total evolution time, in usec. Recommend 10usec.
    """
        
    # Define locations of atoms as a ring
    ground_state_spacing = 1    # How many nearest neighbors the unit disk touches.
                                # 1 is nearest neighbor, 2 is next nearest, etc.
                                # The ground state should be a Z(#+1) state.
    R_interatomic = 9           # The unit disk radius, in um
    L = Nqubits                 # Number of atoms
    
    # Compute the radius and positions of the atoms
    radius = R_interatomic/(2*sin(pi/L * (ground_state_spacing)))
    th = linspace(0,2*pi,L+1)[0:-1]
    coords = array([radius * sin(th) , radius * cos(th)]).T
    
    #print('\nParameters of the atoms:')
    #print('Number of atoms:  {:0.0f}'.format(L))
    #print('Distance of atoms from origin: {:0.6f}um'.format(radius))
    #print('Unit disk radius: {:0.6f}um'.format(R_interatomic))
    
    # Visualize the atoms
    reg = Register.from_coordinates(coords, prefix='atom')
    #reg.draw(blockade_radius=R_interatomic, draw_half_radius=True, draw_graph = True)

    
    # Define the waveform
    
    # The Rydberg energy at the unit disk radius sets the maximum energy scale of delta and omega
    rydberg_blockade_energy = MockDevice.interaction_coeff / R_interatomic**6
    
    Omega_max = 0.9*rydberg_blockade_energy
    
    delta_0 = -2*rydberg_blockade_energy
    delta_f = 0.9*rydberg_blockade_energy
    
    t_rise = 200*Ttotal   # Rise time, in nanoseconds
    t_fall = 200*Ttotal   # Fall time, in nanoseconds
    t_sweep =600*Ttotal   # Sweep time, in nanoseconds
    
    #print('\nParameters of the sweep:')
    #print('Maximum  omega: {:0.6f}'.format(Omega_max))
    #print('Starting delta: {:0.6f}'.format(delta_0))
    #print('Ending   delta: {:0.6f}'.format(delta_f))
    #print('Rise time:  {:0.1f}ns'.format(t_rise))
    #print('Sweep time: {:0.1f}ns'.format(t_sweep))
    #print('Fall time:  {:0.1f}ns'.format(t_fall))
    #print('')
    
    # Define a waveform
    rise = Pulse.ConstantDetuning(RampWaveform(t_rise, 0., Omega_max), delta_0, 0.)
    sweep = Pulse.ConstantAmplitude(Omega_max, RampWaveform(t_sweep, delta_0, delta_f), 0.)
    fall = Pulse.ConstantDetuning(RampWaveform(t_fall, Omega_max, 0.), delta_f, 0.)
    
    # Define a sequence using the pulser interface
    seq = Sequence(reg, MockDevice)
    seq.declare_channel('ising', 'rydberg_global')
    seq.add(rise, 'ising')
    seq.add(sweep, 'ising')
    seq.add(fall, 'ising')
    #seq.draw()
    
    
    
    # Simulate evolution
    #t1 = time.time()
    sim = Simulation(seq, sampling_rate=0.1)
    results = sim.run(progress_bar=False)
    #t2 = time.time()
    #print('Total run time: {:0.4f} seconds'.format(t2-t1))
    
    # Look at the last few elements to compare the state to Bloqade
    print(sort(abs(array(results.get_final_state()))**2,axis=0)[-4::])
    #print('Maximum weight eigenstates should be {:0.0f}-fold degenerate'.format(ground_state_spacing+1))



if __name__=="__main__":
    # Run and time tests
    fil = open('pulser_compare.dat','w')
    nqubits_ = range(4,20)
    for nqubits in nqubits_:
        t1 = time.time()
        dout = pulser_compare(nqubits,10)
        t2 = time.time()
        print('Nqubits: {:2.0f}  |  Total time: {:0.6f}s'.format(nqubits,t2-t1))
        fil.write('{:2.0f}, {:0.14f}\n'.format(nqubits,t2-t1))
        
        if t2-t1>50: # Stop if a time is too long
            break
    fil.close()