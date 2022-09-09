using BloqadeSchema
using Bloqade
using JSON


# task 1. X[pi/2] pulse, wait T, X[pi/2] pulse
# task 2. X[pi/2] pulse, wait T, Y[pi/2] pulse
# task 3. X[pi/2] pulse, wait T/2, X[pi] pulse, wait T/2, X[pi/2] pulse

# task 1

atoms = AtomList([(0.0, 0.0)])

T_pi_2 = 1/16
T= 0.1
epsilon =0.005
Ω = piecewise_constant(clocks=[0.0, epsilon, T_pi_2+ epsilon, T + T_pi_2+epsilon, T + 2* T_pi_2+epsilon, T + 2* T_pi_2+2* epsilon], values= 2π*[0.0, 4, 0, 4, 0])

H = rydberg_h(atoms;Ω=Ω)
h = to_json(H,waveform_tolerance=1e-1,warn=true)


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global5.json","w") do f
    JSON.print(f, h)
end



# task 2

atoms = AtomList([(0.0, 0.0)])

T_pi_2 = 1/16
T= 0.1
epsilon =0.05
Ω = piecewise_constant(clocks=[0.0, epsilon, T_pi_2+ epsilon, T + T_pi_2+epsilon, T + 2* T_pi_2+epsilon, T + 2* T_pi_2+2* epsilon], values= 2π*[0.0, 4, 0, 4, 0])
ϕ = piecewise_constant(clocks=[0.0, epsilon, T_pi_2+ epsilon, T + T_pi_2+epsilon, T + 2* T_pi_2+epsilon, T + 2* T_pi_2+2* epsilon], values= 2π*[0.0, 0, 0, 0, 0])


H = rydberg_h(atoms; Ω=Ω,  ϕ = ϕ)
h = to_json(H,waveform_tolerance=1e-1,warn=true)


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global5.json","w") do f
    JSON.print(f, h)
end

# task 3

atoms = AtomList([(0.0, 0.0)])

T_pi_2 = 1/16
T= 0.1
epsilon =0.005
Ω = piecewise_constant(clocks=[0.0, epsilon, T_pi_2+ epsilon, T/2 + T_pi_2+epsilon, T/2 + 3* T_pi_2+epsilon, T + 3* T_pi_2+epsilon, T + 4* T_pi_2+epsilon, T + 4* T_pi_2+2* epsilon], values= 2π*[0.0, 4, 0, 4, 0, 4, 0])

H = rydberg_h(atoms;Ω=Ω)
h = to_json(H,waveform_tolerance=1e-1,warn=true)


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global5.json","w") do f
    JSON.print(f, h)
end
