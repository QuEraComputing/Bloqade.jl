# Rydberg Blockade

The Rydberg blockade mechanism is a crucial component of the operation of neutral atom computers. By including a term in the effective Hamiltonian which adds an interaction between adjacent Rydberg atoms, we may construct gates and nontrivial dynamics which create entanglement and correlation across the system.

The term is a van der Waals interaction where there is an energy shift if two adjacent atoms are in the Rydberg state. The Van der Waals interaction is

```
V_{ij} = \frac{C_6}{|r_i - r_j|^6}\hat n_i \hat n_j
```

where $\hat n_i=|r_i\rangle\langle r_i|$ is the number operator on the $i$th site, which is 1 if the atom is in the Rydberg state, and zero if the atom is in the ground state. The coefficeint $C_6 = 5393398 MHz/\mu m^6$ is the interaction strength; characteristically, this interaction has a strength $\approx 10MHz$ for two atoms seperated by $10\mu m$

![RydbergBlockade](https://user-images.githubusercontent.com/20091330/156224940-f1558b2e-2762-428d-8cd3-9663d6da3033.png)
