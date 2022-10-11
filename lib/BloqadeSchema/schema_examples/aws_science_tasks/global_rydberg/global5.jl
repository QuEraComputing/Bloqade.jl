using BloqadeSchema
using Bloqade
using JSON

# task 1
# X[pi/2] pulse, wait T, X[pi/2] pulse

atoms = AtomList([(0.0, 0.0)])

T_pi_2 = 1/16
T= 0.1

Ω = piecewise_constant(clocks = [0.0, 
                                 T_pi_2, 
                                 T + T_pi_2, 
                                 T + 2 * T_pi_2], 
                                 values = 2π * [4,0,4])
ϕ = constant(;duration = T + 2 * T_pi_2, value = 0)
Δ = constant(;duration = T + 2 * T_pi_2, value = 0)

h = rydberg_h(atoms;Ω=Ω,ϕ=ϕ,Δ=Δ)
hardware_h, info = hardware_transform(h)
h = to_json(hardware_h) 

open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global5.json","w") do f
    JSON.print(f, JSON.parse(h))
end



# task 2
# X[pi/2] pulse, wait T, Y[pi/2] pulse

atoms = AtomList([(0.0, 0.0)])

T_pi_2 = 1/16
T= 0.1

Ω = piecewise_constant(clocks = [0.0, T_pi_2, T + T_pi_2, T + 2 * T_pi_2], values = 2π * [4,0,4] )
ϕ = piecewise_constant(clocks = [0.0, T_pi_2, T + T_pi_2, T + 2 * T_pi_2], values =  [0.0,0.0,π/2] )
Δ = constant(;duration = T + 2 * T_pi_2, value = 0)

h = rydberg_h(atoms;Ω=Ω,ϕ=ϕ,Δ=Δ)
hardware_h, info = hardware_transform(h)
h = to_json(hardware_h) 


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global5.json","w") do f
    JSON.print(f, JSON.parse(h))
end

# task 3
# X[pi/2] pulse, wait T/2, X[pi] pulse, wait T/2, X[pi/2] pulse

atoms = AtomList([(0.0, 0.0)])

T_pi_2 = 1/16
T= 0.1

Ω = piecewise_constant(clocks=[0.0, 
                               T_pi_2,
                               T/2 + T_pi_2, 
                               T/2 + 3*T_pi_2, 
                               T + 3*T_pi_2,
                               T + 4*T_pi_2], 
                               values = 2π * [4, 0, 4, 0, 4])
ϕ = constant(;duration=T+4*T_pi_2, value = 0)
Δ = constant(;duration=T+4*T_pi_2, value = 0)

h = rydberg_h(atoms;Ω=Ω,ϕ=ϕ,Δ=Δ)
h_hardware, info = hardware_transform(h)
h = to_json(h_hardware)


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global5.json","w") do f
    JSON.print(f, JSON.parse(h))
end
