using OrdinaryDiffEq
struct ODEState{VT}
    v::VT
    t::Float64
end
Base.copy(s::ODEState) = ODEState(copy(s.v), s.t)

function evolve_step!(state::ODEState; integrator, dt)
    step!(integrator, dt)
    return ODEState(state.v, state.t + dt)
end

function field_function!(du, u, p, t)
    du .= sin(p[1]*t)
    return du
end
u0 = zeros(3)
tspan = (0.0, 2.0)

prob = ODEProblem(field_function!, u0, tspan, [2.0])
integrator = init(prob, Vern8(); dt=0.01)
step!(integrator, 0.01)

# t
# dt
# uprev
# destats.nf  # number of fevals

# fixed
# ---------------------
# alg    # algorithm
# opts   # options
# f
# 
# not useful
# ---------------------
# u  # overwritten
# p  # parameters
# EEst   # error estimation, not used.
# eigen_est # ?, not used

# k      # f values.

# cache.tab.extra
# k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,tmp,rtmp,atmp

# option ref: https://diffeq.sciml.ai/stable/basics/common_solver_opts/#solver_options
function diffode(; alg=Vern8(),
        dt=0.01, alg_hint=:nonstiff)
    prob = ODEProblem(f, u0, tspan)
    integrator = init(prob, alg; dt, alg_hint)
    step!(integrator, dt)
end