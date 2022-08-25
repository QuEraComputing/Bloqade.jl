# measurements.jl
#
# Defines estimators and provides measurements

function sample(H::AbstractIsing, qmc_state::BinaryQMCState, M::Int=length(qmc_state.operator_list) ÷ 2)
    operator_list = qmc_state.operator_list
    spin_prop = copy!(qmc_state.propagated_config, qmc_state.left_config)

    @inbounds for i in 1:M
        op = operator_list[i]
        if !isdiagonal(H, op)
            spin_prop[getsite(H, op)] ⊻= 1 #spinflip
        end
    end
    return spin_prop
end

function update_two_pt_fn!(two_pt_fn::Matrix{T}, one_pt_fn::Vector{T}, sample::Vector, n::Int) where T <: AbstractFloat
    n += 1
    delta = @. T(sample) - one_pt_fn

    LinearAlgebra.BLAS.syr!('U', inv(n), delta, two_pt_fn)  # two_pt_fn .+= (delta/n) * transpose(delta)
    two_pt_fn *= (n-1)/n

    one_pt_fn .+= delta / n

    return two_pt_fn, one_pt_fn, n
end

function simulation_cell(H::AbstractIsing, qmc_state::BinaryQMCState, r::OrdinalRange{Int, Int})
    operator_list = qmc_state.operator_list

    cell = falses(length(qmc_state.left_config), length(r))
    spin_prop = copy(qmc_state.left_config)
    c = 1

    @inbounds for (n, op) in enumerate(operator_list)
        if !isdiagonal(H, op)
            spin_prop[getsite(H, op)] ⊻= 1 #spinflip
        end
        if n in r
            copy!(view(cell, :, c), spin_prop)
            c += 1
        end
    end
    return cell
end
simulation_cell(H::AbstractIsing, qmc_state::BinaryQMCState) = simulation_cell(H, qmc_state, 1:length(qmc_state.operator_list))

# 1 -> spin-up (+1)
# 0 -> spin-down (-1)
magnetization(spin_prop) = 2*mean(spin_prop) - 1

function staggered_magnetization(H::AbstractRydberg, spin_prop)
    M = 0.0

    if H.lattice isa Rectangle
        spin_prop = reshape(spin_prop, H.lattice.n1, H.lattice.n2)

        for j in axes(spin_prop, 2), i in axes(spin_prop, 1)
            M += ((-1)^(i + j)) * (spin_prop[i, j] - 0.5)
        end
    else
        for i in eachindex(spin_prop)
            M += ((-1) ^ i) * (spin_prop[i] - 0.5)
        end
    end

    return M / nspins(H)
end

function kagome_nematic(lattice::Kagome, sublattice::Vector{Int}, spin_prop)
    # https://www.pnas.org/content/pnas/118/4/e2015785118.full.pdf
    # page 4

    N2 = lattice.n2 * 2
    Nc = 3 * (N2^2)
    phase = exp((2*π/3)*im)
    spin_prop_A = spin_prop[findall(x -> x == 1, sublattice)]
    spin_prop_B = spin_prop[findall(x -> x == 2, sublattice)]
    spin_prop_C = spin_prop[findall(x -> x == 0, sublattice)]

    Φ = sum(spin_prop_A) + phase * sum(spin_prop_B) + phase^2 * sum(spin_prop_C)
    Φ *= 3/Nc

    return Φ
end

function correlation_functions(spin_prop)
    N = size(spin_prop)[1]

    correlations = zeros(N, N)
    for i in 1:(N-1)
        for j in (i+1):N
            correlations[i,j] = spin_prop[i] * spin_prop[j]
        end
    end

    return correlations
end


function domain_wall_density(H::AbstractRydberg, spin_prop)
    # Currently set up for 1D + open boundary conditions
    #=
    From https://www.nature.com/articles/nature24622:

    Domain walls are identified as either two neighbouring atoms
    in the same state or a ground-state atom at the edge of the array.
    =#

    L = nspins(H)
    dwd = 0.0

    # check boundaries
    dwd += iszero(spin_prop[1]) ? 1.0 : 0.0
    dwd += iszero(spin_prop[end]) ? 1.0 : 0.0

    # now check the bulk
    for i in 1:(L-1)
        dwd += (spin_prop[i] == spin_prop[i+1]) ? 1.0 : 0.0
    end

    return dwd / L
end

num_single_site_diag(H::AbstractIsing, operator_list) = mean(x -> issiteoperator(H, x) && isdiagonal(H, x), operator_list)
num_single_site_offdiag(H::AbstractIsing, operator_list) = mean(x -> issiteoperator(H, x) && !isdiagonal(H, x), operator_list)
num_single_site(H::AbstractIsing, operator_list) = mean(issiteoperator(H), operator_list)
num_two_site_diag(H::AbstractIsing, operator_list) = mean(isbondoperator(H), operator_list)


function autocorrelation(m::Vector)
    N = length(m)

    m′ = m .- mean(m)
    m′ = vcat(m′, zeros(N))
    mw = fft(m′)
    s = abs2.(mw)

    @inbounds chi = real(ifft(s)[1:N])

    @inbounds for i in 1:N
        chi[i] /= (2*N)  # normalize FFT
        chi[i] /= (N - i - 1)
    end
    return chi
end


# use method explained by Sokal to estimate correlation time
# https://pdfs.semanticscholar.org/0bfe/9e3db30605fe2d4d26e1a288a5e2997e7225.pdf
function correlation_time(m::Vector)
    ac = autocorrelation(m)
    ac_0 = ac[1]

    corr_time = 0.0
    @inbounds for M in axes(ac, 1)
        corr_time += (ac[M] / ac_0)
        if M >= 10*corr_time
            break
        end
    end

    return corr_time
end
