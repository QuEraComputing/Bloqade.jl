export cmaes
include("CMAESCore.jl")

# Acknowledgement:
#
# This code is adapted from: https://github.com/wildart/Evolutionary.jl
#
# Covariance Matrix Adaptation Evolution Strategy
# ==================
#
# Implementation: (μ/μ_I,λ)-CMA-ES
#
# μ is the number of parents
# λ is the number of offspring.
#

############################ New Interface ############################
# * less input parameters, like number of iterations, tol, verbose.
# * clean `call_back` support, with full access to runtime information.

struct GeneralOptProblem{FT}
	lossfunc::FT
	N::Int
end
ellength(opt::GeneralOptProblem) = opt.N
loss(opt::GeneralOptProblem, x) = opt.lossfunc(x)

struct CMAESIter{PT, T}
	cr::CMAESRuntime{T}
	prob::PT
end
function CMAESIter(prob, individual; num_parents::Integer, num_offsprings::Integer, σ::Real=1)
	cr = CMAESRuntime(individual, num_parents=num_parents, num_offsprings=num_offsprings, N=ellength(prob), σ=σ)
	CMAESIter(cr, prob)
end

function CMAESIter(lossfunc::Function, individual; num_parents::Integer, num_offsprings::Integer, σ::Real=1)
	prob = GeneralOptProblem(lossfunc, length(individual))
	CMAESIter(prob, individual, num_parents=num_parents, num_offsprings=num_offsprings, σ=σ)
end

function Base.iterate(ci::CMAESIter, state=1)
	cmaes_step!(ci.cr, ci.prob, τ=init_τ(ci.prob), τ_c=init_τ_c(ci.prob), τ_σ=init_τ_σ(ci.prob)), state+1
end

####################### Problem definition ###############################
function populate!(fitoff::AbstractVector, offspring::AbstractVector, prob)
    for i in 1:length(fitoff)
        fitoff[i] = loss(prob, offspring[i]) # Evaluate fitness
    end
end


"""
    cmaes(f, x0; num_offsprings=50, num_parents=20, maxiter=1000, tol=1e-8)

`CMA-ES` optimizer.
"""
function cmaes(f, x0; num_offsprings=50, num_parents=20, maxiter=1000, tol=1e-8)
	ci = CMAESIter(f, x0, num_offsprings=num_offsprings, num_parents=num_parents)
	for (count, cr) in enumerate(ci)
		if count == maxiter || cr.σ<tol break end
		println("BEST: $(cr.fitpop[1]): $(cr.σ)")
	end
	result, fitness = best(ci.cr)
    return result
end
