using LinearAlgebra

function evolution_operator(t::Float64)
    ϕ = 2.2 * sin(π * t)^2
    U = zeros(ComplexF64, 2,2)
    U[1,1] =  1 / sqrt(2)
    U[2,1] =  1 / sqrt(2)
    U[2,2] =  1 / sqrt(2)
    U[1,2] = -1 / sqrt(2)

    U * diagm(exp.([-im*ϕ, im*ϕ])) * U'
end

function solution(t)
    U = evolution_operator(t)
    return U * [1.0, 0.0]

end
