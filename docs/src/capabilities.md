# [Hardware Capabilities](@id capabilities)

Although tools such as those that Bloqade offers in its [Schema](@ref schema) can automatically transform user-created Hamiltonians to conform to [QuEra's Quantum Processor *Aquila*](https://www.quera.com/aquila) requirements with minimal user intervention, the following information may still be helpful in designing Hamiltonians to run on hardware.

## Programmatic Access

You can obtain the values from the human-readable form below programmatically through the [`BloqadeSchema.get_device_capabilities`](@ref) function which returns `BloqadeSchema.DeviceCapabilities`, a named tuple.

```@repl
using BloqadeSchema
device_capabilities = get_device_capabilities();
BloqadeSchema.pprint(device_capabilities)
device_capabilities.rydberg.global_value.detuning_max # get the maximum global detuning value
```

The field names for each value have been provided below but will require you to provide the proper prefix like in the example above (e.g. the Minimum Rabi Frequency is under `rabi_frequency_min` which exists under `global_value` under `rydberg` in the `DeviceCapabilities` tuple, and should be access via *`your_device_capabilities_instance`*`.rydberg.global_value.rabi_frequency_min`)

## *Aquila* Capabilities

### Task

| Capability              | Field              | Value |
|:------------------------|:-------------------|:------|
| Minimum Number of Shots | `number_shots_min` | 1     |
| Maximum Number of Shots | `number_shots_max` | 1000  |

### Lattice Geometry

| Capability                              | Field                  | Value   |
|:----------------------------------------|:-----------------------|:--------|
| Maximum Number of Qubits                | `number_qubits_max`    | 256     |
| Maximum Lattice Area Width              | `width`                | 75.0 µm |
| Maximum Lattice Area Height             | `height`               | 76.0 µm |
| Minimum Radial Spacing between Qubits   | `spacing_radial_min`   | 4.0 µm  |
| Minimum Vertical Spacing between Qubits | `spacing_vertical_min` | 4.0 µm  |
| Position Resolution                     | `position_resolution`  | 0.1 µm  |
| Maximum Number of Sites                 | `number_sites_max`     | 256     |

### Global Rydberg Values

| Capability                       | Field                          | Value              |
|:---------------------------------|:-------------------------------|:-------------------|
| Rydberg Interaction Constant     | `c6_coefficient`               | 5.42×10⁶ MHz × µm⁶ |
| Minimum Rabi Frequency           | `rabi_frequency_min`           | 0.00 MHz           |
| Maximum Rabi Frequency           | `rabi_frequency_max`           | 15.8 MHz           |
| Rabi Frequency Resolution        | `rabi_frequency_resolution`    | 0.0004 MHz         |
| Maximum Rabi Frequency Slew Rate | `rabi_frequency_slew_rate_max` | 250.0 MHz/µs       |
| Minimum Detuning                 | `detuning_min`                 | -125.0 MHz         |
| Maximum Detuning                 | `detuning_max`                 | 125.0 MHz          |
| Detuning Resolution              | `detuning_resolution`          | 2.0×10⁻⁷ MHz       |
| Maximum Detuning Slew Rate       | `detuning_slew_rate_max`       | 2500.0 MHz/µs      |
| Minimum Phase                    | `phase_min`                    | -99.0 rad          |
| Maximum Phase                    | `phase_max`                    | 99.0 rad           |
| Phase Resolution                 | `phase_resolution`             | 5.0×10⁻⁷ rad       |
| Minimum Time                     | `time_min`                     | 0.0 µs             |
| Maximum Time                     | `time_max`                     | 4.0 µs             |
| Time Resolution                  | `time_resolution`              | 0.001 µs           |
| Minimum Δt                       | `time_delta_min`               | 0.05 µs            |