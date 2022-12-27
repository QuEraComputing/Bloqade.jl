# [Hardware Capabilities](@id capabilities)

Although tools such as those that Bloqade offers in its [Schema](@ref schema) can automatically transform user-created Hamiltonians to conform to [QuEra's Quantum Processor *Aquila*](https://www.quera.com/aquila) requirements with minimal user intervention, the following information may still be helpful in designing Hamiltonians to run on hardware.

### Unit Conversions

The values below have been converted from their base SI unit representations (viewable in the `lib/BloqadeSchema/config/capabilities-qpu1-mock.json` file of the repository with units specified by `lib/BloqadeSchema/config/capabilities-qpu1-mock-units.json`) to the ones Bloqade uses by default. To display the converted results, the values have been rounded to two decimal places accompanied by the "≈" sign.

## *Aquila* Capabilities

### Task

| Capability              | Value |
|:------------------------|:------|
| Minimum Number of Shots | 1     |
| Maximum Number of Shots | 1000  |

### Lattice Geometry

| Capability                              | Value  |
|:----------------------------------------|:-------|
| Maxmimum Number of Atoms                | 256    |
| Maximum Lattice Area Width              | 75 µm  |
| Maximum Lattice Area Height             | 76 µm  |
| Minimum Radial Spacing between Qubits   | 4.0 µm |
| Minimum Vertical Spacing between Qubits | 4.0 µm |
| Position Resolution                     | 0.1 µm |
| Maximum Number of Sites                 | 256    |

### Global Rydberg Values

| Capability                       | Value                      |
|:---------------------------------|:---------------------------|
| Rydberg Interaction Constant     | ≈8.63×10⁻³¹ 2π × MHz × µm⁶ |
| Minimum Rabi Frequency           | 0.00 2π × MHz              |
| Maximum Rabi Frequency           | ≈2.51 2π × MHz             |
| Rabi Frequency Resolution        | ≈6.37×10⁻⁵ 2π × MHz        |
| Maximum Rabi Frequency Slew Rate | ≈39.79 2π × MHz/µs         |
| Minimum Detuning                 | ≈-19.89 2π × MHz           |
| Maximum Detuning                 | ≈19.89 2π × MHz            |
| Detuning Resolution              | ≈3.18×10⁻⁸ 2π × MHz        |
| Maximum Detuning Slew Rate       | ≈397.89 2π × MHz/µs        |
| Minimum Phase                    | 99.0 rad                   |
| Maximum Phase                    | 99.0 rad                   |
| Phase Resolution                 | 5.0×10⁻⁷ rad               |
| Minimum Time                     | 0.0 µs                     |
| Maximum Time                     | 4.0 µs                     |
| Time Resolution                  | 0.001 µs                   |
| Minimum Time Delta               | 0.05 µs                    |