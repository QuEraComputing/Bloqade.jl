# [Hardware Capabilities](@id capabilities)

Although tools such as those that Bloqade offers in its [Schema](@ref schema) can automatically transform user-created Hamiltonians to conform to [QuEra's Quantum Processor *Aquila*](https://www.quera.com/aquila) requirements with minimal user intervention, the following information may still be helpful in designing Hamiltonians to run on hardware.

## *Aquila* Capabilities

### Task

| Capability              | Value |
|:------------------------|:------|
| Minimum Number of Shots | 1     |
| Maximum Number of Shots | 1000  |

### Lattice Geometry

| Capability                              | Value   |
|:----------------------------------------|:--------|
| Maximum Number of Atoms                 | 256     |
| Maximum Lattice Area Width              | 75.0 µm |
| Maximum Lattice Area Height             | 76.0 µm |
| Minimum Radial Spacing between Qubits   | 4.0 µm  |
| Minimum Vertical Spacing between Qubits | 4.0 µm  |
| Position Resolution                     | 0.1 µm  |
| Maximum Number of Sites                 | 256     |

### Global Rydberg Values

| Capability                       | Value              |
|:---------------------------------|:-------------------|
| Rydberg Interaction Constant     | 5.42×10⁶ MHz × µm⁶ |
| Minimum Rabi Frequency           | 0.00 MHz           |
| Maximum Rabi Frequency           | 15.8 MHz           |
| Rabi Frequency Resolution        | 0.0004 MHz         |
| Maximum Rabi Frequency Slew Rate | 250.0 MHz/µs       |
| Minimum Detuning                 | -125.0 MHz         |
| Maximum Detuning                 | 125.0 MHz          |
| Detuning Resolution              | 2.0×10⁻⁷ MHz       |
| Maximum Detuning Slew Rate       | 2500.0 MHz/µs      |
| Minimum Phase                    | -99.0 rad          |
| Maximum Phase                    | 99.0 rad           |
| Phase Resolution                 | 5.0×10⁻⁷ rad       |
| Minimum Time                     | 0.0 µs             |
| Maximum Time                     | 4.0 µs             |
| Time Resolution                  | 0.001 µs           |
| Minimum Δt                       | 0.05 µs            |