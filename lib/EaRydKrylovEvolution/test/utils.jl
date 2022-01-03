using Graphs
using EaRydKrylovEvolution

const test_graph = SimpleGraph(5)
add_edge!(test_graph, 1, 2)
add_edge!(test_graph, 2, 3)
add_edge!(test_graph, 2, 4)
add_edge!(test_graph, 2, 5)
add_edge!(test_graph, 3, 4)
add_edge!(test_graph, 4, 5)

const test_subspace_v = [0, 1, 2, 4, 5, 8, 9, 16, 17, 20, 21]
const test_subspace = independent_set_subspace(test_graph)

create_test_hamiltonian(Δ, Ω, ϕ) = create_test_hamiltonian!(zeros(ComplexF64, 11, 11), Δ, Ω, ϕ)

function create_test_hamiltonian!(H::AbstractMatrix, Δ, Ω, ϕ)
    Ω = Ω ./ 2
    H[1,1] = Δ[1] + Δ[2] + Δ[3] + Δ[4] + Δ[5]
    H[2,2] = - Δ[1] + Δ[2] + Δ[3] + Δ[4] + Δ[5]
    H[3,3] = Δ[1] - Δ[2] + Δ[3] + Δ[4] + Δ[5]
    H[4,4] = Δ[1] + Δ[2] - Δ[3] + Δ[4] + Δ[5]
    H[5,5] = - Δ[1] + Δ[2] - Δ[3] + Δ[4] + Δ[5]
    H[6,6] = Δ[1] + Δ[2] + Δ[3] - Δ[4] + Δ[5]
    H[7,7] = - Δ[1] + Δ[2] + Δ[3] - Δ[4] + Δ[5]
    H[8,8] = + Δ[1] + Δ[2] + Δ[3] + Δ[4] - Δ[5]
    H[9,9] = - Δ[1] + Δ[2] + Δ[3] + Δ[4] - Δ[5]
    H[10,10] = + Δ[1] + Δ[2] - Δ[3] + Δ[4] - Δ[5]
    H[11,11] = - Δ[1] + Δ[2] - Δ[3] + Δ[4] - Δ[5]

    # spin flip terms
    H[1,2] = Ω[1] * exp(im * ϕ[1])
    H[2,1] = Ω[1] * exp(-im * ϕ[1])

    H[1,3] = Ω[2] * exp(im * ϕ[2])
    H[3,1] = Ω[2] * exp(-im * ϕ[2])

    H[1,4] = Ω[3] * exp(im * ϕ[3])
    H[4,1] = Ω[3] * exp(-im * ϕ[3])

    H[1,6] = Ω[4] * exp(im * ϕ[4])
    H[6,1] = Ω[4] * exp(-im * ϕ[4])

    H[1,8] = Ω[5] * exp(im * ϕ[5])
    H[8,1] = Ω[5] * exp(-im * ϕ[5])

    H[2,5] = Ω[3] * exp(im * ϕ[3])
    H[5,2] = Ω[3] * exp(-im * ϕ[3])

    H[2,7] = Ω[4] * exp(im * ϕ[4])
    H[7,2] = Ω[4] * exp(-im * ϕ[4])

    H[2,9] = Ω[5] * exp(im * ϕ[5])
    H[9,2] = Ω[5] * exp(-im * ϕ[5])

    H[4,5] = Ω[1] * exp(im * ϕ[1])
    H[5,4] = Ω[1] * exp(-im * ϕ[1])

    H[4,10] = Ω[5] * exp(im * ϕ[5])
    H[10,4] = Ω[5] * exp(-im * ϕ[5])

    H[5,11] = Ω[5] * exp(im * ϕ[5])
    H[11,5] = Ω[5] * exp(-im * ϕ[5])

    H[6,7] = Ω[1] * exp(im * ϕ[1])
    H[7,6] = Ω[1] * exp(-im * ϕ[1])

    H[8,9] = Ω[1] * exp(im * ϕ[1])
    H[9,8] = Ω[1] * exp(-im * ϕ[1])

    H[8,10] = Ω[3] * exp(im * ϕ[3])
    H[10,8] = Ω[3] * exp(-im * ϕ[3])

    H[8,10] = Ω[3] * exp(im * ϕ[3])
    H[10,8] = Ω[3] * exp(-im * ϕ[3])

    H[9,11] = Ω[3] * exp(im * ϕ[3])
    H[11,9] = Ω[3] * exp(-im * ϕ[3])

    H[10,11] = Ω[1] * exp(im * ϕ[1])
    H[11,10] = Ω[1] * exp(-im * ϕ[1])
    return H
end
