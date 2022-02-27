using OrdinaryDiffEq
struct ODEState{VT}
    u::VT
    t::Float64
end
Base.copy(s::ODEState) = ODEState(copy(s.u), s.t)

function evolve_step!(state::ODEState; integrator, dt)
    step!(integrator, dt)
    return ODEState(state.v, state.t + dt)
end

function field_function!(du, u, p, t)
    du .= sin(p[1]*t)
    return du
end
function save_state(integrator)
    return ODEState(copy(integrator.u), integrator.t)
end

function load_state!(integrator, state)
    integrator.u .= state.u
    integrator.t = state.t
    return integrator
end

function test_save_load()
    integrator = init(prob, Vern8(); dt=0.01)
    step!(integrator, 0.01, true)
    step!(integrator, 0.01, true)
    d = save_state(integrator)
    step!(integrator, 0.01, true)
    step!(integrator, 0.01, true)
    u1 = copy(integrator.u)
    load_state!(integrator, d)
    step!(integrator, 0.01, true)
    step!(integrator, 0.01, true)
    u2 = copy(integrator.u)
    u1 ≈ u2
end

using TreeverseAlgorithm
function integrator_step(integrator, dt)
    function (state)
        step!(load_state!(integrator, state) , dt, true)
        return save_state(integrator)
    end
end

function integrator_back(integrator, dt)
    function (sin, gout)
        if gout === nothing
            @show sin
            return fill(1.0, 3), 0.0, 0.0
        else
            gsout, gt, gp = gout
            p = integrator.p[1]
            return (gsout, sum(gsout .* (p * cos(sin.t * p))), gp .+ sum(gsout .* dt * sin.t * cos(sin.t*p)))
        end
    end
end

dt = 0.01
u0 = zeros(3)
tspan = (0.0, 2.0)
prob = ODEProblem(field_function!, u0, tspan, [2.0])
integrator = init(prob, Vern8(); dt=dt)
treeverse(integrator_step(integrator, dt), integrator_back(integrator, dt), ODEState(u0, 0.0); f_inplace=false, N=200, δ=5)

# t
# dt
# uprev
# k      # f values.

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

# destats.nf  # number of fevals

# cache.tab.extra
# k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,tmp,rtmp,atmp

# option ref: https://diffeq.sciml.ai/stable/basics/common_solver_opts/#solver_options
function diffode(; alg=Vern8(),
        dt=0.01, alg_hint=:nonstiff)
    prob = ODEProblem(f, u0, tspan)
    integrator = init(prob, alg; dt, alg_hint)
    step!(integrator, dt)
end