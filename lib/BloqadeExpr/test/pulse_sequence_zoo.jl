using BloqadeExpr
using BloqadeODE, BloqadeKrylov
using Yao
using BitBasis
using Test

function single_site_rydberg_params(α, θ, ϕ; t = pi)
    Δ = cos(θ) * α / t
    Ω = sqrt((α/t)^2 - Δ^2)
    ϕ_0 = -ϕ
    return Ω, ϕ_0, Δ
end

function single_site_pulse_params(Ω, ϕ_0, Δ; t = pi)
    α = sqrt(abs2(Ω) + Δ^2) * t
    θ = acos(Δ/α*t)
    ϕ = -ϕ_0
    return (α, θ, ϕ)
end

function single_site_hamiltonian(atoms, j, α, θ, ϕ; name = :rydberg, t = π)
    Ω, ϕ, Δ = single_site_rydberg_params(α, θ, ϕ; t)
    n = length(atoms)
    @assert j in 1:n
    mask = zeros(n)
    mask[j] = 1
    Ω = Ω * mask
    ϕ = ϕ * mask
    Δ = Δ * mask
    name === :rydberg && (h = rydberg_h_3(atoms; Ω_r = Ω, ϕ_r = ϕ, Δ_r = Δ))
    name === :hyperfine && (h = rydberg_h_3(atoms; Ω_hf = Ω, ϕ_hf = ϕ, Δ_hf = Δ))
    return h
end

function single_site_evolution(reg, atoms, j, α, θ, ϕ; name = :rydberg, t = π)
    h = single_site_hamiltonian(atoms, j, α, θ, ϕ; name, t)
    return KrylovEvolution(reg, 0.0:1e-3:t, h)
end

@testset "11-pulse SWAP" begin
    # Notes that this SWAP only acts correctly on |0⟩ and |1⟩. For example, |01⟩ → |10⟩ but |0r⟩ → |1r⟩.
    atoms = [(0.0, 0.0), (0.0, 1)]
    sequence = [(:rydberg, 2, [π, π / 2, 0]),
        (:hyperfine, 2, [π, π / 2, 0]),
        (:rydberg, 1, [π, π / 2, 0]),
        (:hyperfine, 1, [π, π / 2, 0]),
        (:rydberg, 2, [-π, π / 2, 0]),
        (:rydberg, 1, [-π, π / 2, 0]),
        (:hyperfine, 1, [-π, π / 2, 0]),
        (:rydberg, 2, [π, π / 2, 0]),
        (:hyperfine, 2, [-π, π / 2, 0]),
        (:rydberg, 1, [π, π / 2, 0]),
        (:rydberg, 2, [-π, π / 2, 0])]
    hs = [single_site_hamiltonian(atoms, j, ps...; name) for (name, j, ps) in sequence]
    ms = [Matrix(h) for h in hs]
    @test length(findall(x->!isapprox(x, 0;atol = 1e-6), prod(exp(m*im*pi) for m in ms))) == 9
    ids = [1, 2, 4, 5]
    SW = matblock(prod(exp(-m*im*pi) for m in ms)[ids, ids])
    @test operator_fidelity(SW, SWAP) > 1 - 1e-8
    reg = product_state(dit"10;3")
    evos = [single_site_evolution(reg, atoms, j, ps...; name) for (name, j, ps) in sequence]
    for evo in evos
        emulate!(evo)
    end
    @test isapprox(abs(state(reg)[2]), 1; atol = 1e-6)
end

function global_hamiltonian(atoms, ξ, Δ; name = :rydberg)
    name === :rydberg && (h = rydberg_h_3(atoms; Ω_r = 1.0, ϕ_r = ξ, Δ_r = Δ))
    name === :hyperfine && (h = rydberg_h_3(atoms; Ω_hf = 1.0, ϕ_hf = ξ, Δ_hf = Δ))
    return h
end

function global_krylov_evolution(reg, atoms, ξ, Δ, τ; name = :rydberg, step = 1e-3)
    h = global_hamiltonian(atoms, ξ, Δ; name)
    return KrylovEvolution(reg, 0:step:τ, h)
end

function global_ode_evolution(reg, atoms, ξ, Δ, τ; name = :rydberg)
    h = global_hamiltonian(atoms, ξ, Δ; name)
    return SchrodingerProblem(reg, (0.0, τ), h)
end

function test_global_pulse(reg, atoms, sequence; op = SWAP, atol = 1e-6)
    goal = copy(state(reg))
    M = nothing
    for (name, ps) in sequence
        h = global_hamiltonian(atoms, ps[1:2]...; name)
        t = ps[3]
        mh = exp(-im*t*Matrix(h))
        isnothing(M) ? (M = mh) : (M = mh*M)
        goal = mh * goal
        evo = SchrodingerProblem(reg, (0, t), h)
        emulate!(evo)
        fdlt = abs((goal' * state(reg))[1,1])
        @test isapprox(fdlt, 1; atol)
    end
    if length(atoms) == 2
        ids = [1,2,4,5]
        @test operator_fidelity(matblock(M[ids, ids]), op) > 1 - atol
    end
end

@testset "9-pulse SWAP" begin
    atoms = [(0.0, 0.0), (1.0, 0.0)]
    Δ = 0.377371
    ξ = 3.90242
    τ = 4.29268
    xξ = 0.42072
    xΔ = 0.37131
    xτ = 1.60222
    sequence = [
        (:hyperfine, [xξ, xΔ, xτ]),
        (:rydberg, [0.0, Δ, τ]),
        (:rydberg, [ξ, Δ, τ]),
    ]
    sequence = [sequence; sequence; sequence]
    
    atoms = [(0.0, 0.0), (4.0, 0.0)]
    reg = product_state(dit"01;3")
    test_global_pulse(reg, atoms, sequence; op = SWAP, atol = 1e-6)
end

# end

@testset "8-pulse SWAP" begin
    sequence = [
        (:hyperfine, [0.0, -0.39308582446292734, 5.856543874719311]),
        (:rydberg, [-1.373301896192175, 0.28072928877143755, 6.80710127878679]),
        (:rydberg, [1.8110512905361198, 0.06474366573998236, 3.8636876491429755]),
        (:hyperfine, [-4.4836982638809975, -0.0003415966053856353, -3.1341957440477866]),
        (:rydberg, [3.5106267226727317, 0.8048234980426018, 3.8704640030518624]),
        (:rydberg, [1.2538197255706265, 0.8281666963495008, 3.8449644067823194]),
        (:hyperfine, [1.3687393944450619, 0.0056118212964182, 3.1321379883863814]),
        (:rydberg, [-2.7362167093986796, -0.0027837111779557364, -3.1331303677146596])
    ]
    
    atoms = [(0.0, 0.0), (4.0, 0.0)]
    reg = product_state(dit"01;3")
    test_global_pulse(reg, atoms, sequence; op = SWAP, atol = 1e-3)
end

@testset "Inexact 7-pulse SWAP" begin
    sequence = [
        (:rydberg, [0.0, 0.022052631957254598, 9.175649911608547]),
        (:hyperfine, [-1.3068549602703559, -0.0019116955457488213, 3.125447757701348]),
        (:rydberg, [-3.7701937764746623, 0.0004003811953280955, 4.442387518897087]),
        (:hyperfine, [2.400509785176155, 0.003397804044057438, 3.1294709767549667]),
        (:rydberg, [4.543440853588044, -0.5476584731820732, 4.349235209223257]),
        (:rydberg, [4.274590180638748, 0.0740621269411009, 4.250631386650806]),
        (:hyperfine, [1.1589994986561472, -0.33664894310734195, -5.954943545548953])
    ]
    
    atoms = [(0.0, 0.0), (4.0, 0.0)]
    reg = product_state(dit"01;3")
    test_global_pulse(reg, atoms, sequence; op = SWAP, atol = 1e-2)
end

@testset "Exact 7-pulse SWAP" begin
    sequence = [
        (:rydberg, [0.0, -0.1813167182713958, -3.6968102973773536]),
        (:rydberg, [-0.5417186274378146, 2.348000206150177, -1.608814643831605]),
        (:hyperfine, [0.0, 0.0, pi]),
        (:rydberg, [1.657341528146135, -0.8164965718969145, -3.847649512868359]),
        (:rydberg, [-0.5795339591221397, -0.8164965718969145, -3.847649512868359]),
        (:hyperfine, [2.5349178162096853, 0.0, pi]),
        (:rydberg, [3.0349276621572367, 0.0, pi])
    ]
    
    atoms = [(0.0, 0.0), (4.0, 0.0)]
    reg = product_state(dit"01;3")
    test_global_pulse(reg, atoms, sequence; op = SWAP, atol = 1e-6)
end

@testset "Analytic 7-pulse SWAP" begin
    Δ = 0.3644747433308944
    τ = 4.302298535095522
    ξ = -5.537047938127648
    ξ₅ = -π / 4
    ξ₆ = 3π / 4
    sequence = [
        (:rydberg, [0.0, Δ, τ]),
        (:rydberg, [ξ, Δ, τ]),
        (:hyperfine, [0.0, 0.0, π]),
        (:rydberg, [0.0, 0.0, 2π / √2]),
        (:hyperfine, [ξ₅, 0.0, π]),
        (:rydberg, [ξ₆ + ξ, Δ, -τ]),
        (:rydberg, [ξ₆, Δ, -τ]),
    ]

    atoms = [(0.0, 0.0), (4.0, 0.0)]
    reg = product_state(dit"01;3")
    test_global_pulse(reg, atoms, sequence; op = SWAP, atol = 1e-6)
end

@testset "3-pulse CZ" begin
    sequence = [
        (:rydberg, [0.0, 0.8126945875767624, 2.0613620704810827]),
        (:rydberg, [6.153382130563536, -0.4640923475669967, 3.9000469759842877]),
        (:rydberg, [6.277208126588368, 0.8132262498823576, 1.8627674483306587])
    ]
    op = matblock(Matrix(Diagonal([0.9999999999997726 + 0.0im; 
        -0.4306510538860446 + 0.9025185149273534im; 
        -0.4306510538860447 + 0.9025185149273535im; 
        0.6278987469593043 + 0.7782949718113085im])))
    atoms = [(0.0, 0.0), (4.0, 0.0)]
    reg = product_state(dit"11;3")
    test_global_pulse(reg, atoms, sequence; op, atol = 1e-7)
end

@testset "3-pulse CCZ" begin
    sequence = [
        (:rydberg, [0.0, -0.13279501943121658, 6.45885351211801]),
        (:rydberg, [0.7133166536257461, 0.9783159812194855, 4.08145583869183]),
        (:rydberg, [1.4266329023945055, -0.13279514815653518, 6.458853255668704])
    ]
    atoms = [(0.0, 0.0), (2*sqrt(3), 2), (-2*sqrt(3), 2)]
    reg = product_state(dit"111;3")
    test_global_pulse(reg, atoms, sequence; op = control(3, (2, 3), 1=>Z), atol = 1e-9)
    state(reg)
end