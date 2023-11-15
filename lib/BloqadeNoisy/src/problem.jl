"""
    struct NoisySchrodingerEquation

Decorator object for SchrodingerEquation. Implements the f(dstate, state, p, t)
interface to fit into a standard ODE Problem.

# Arguments

- `equation`: The coherent 'SchrodingerEquation` equation for the problem
- `imag_evo`: ∑_L L^†L for collapse operators L which defines the non-hermitian part of the Hamiltonian
- `num_cops``: number of collapse operators (for displaying)
"""
struct NoisySchrodingerEquation{mType}
    equation::SchrodingerEquation #Wraps this SchrodingerEquation
    imag_evo::mType #anti-hermitian part of the Hamiltonian
    num_cops::Int
end

function (eq::NoisySchrodingerEquation)(dstate, state, p, t::Number)
    eq.equation(dstate, state, p, t) #d/dt|ψ> = -iH|ψ> - .5∑_L L^†L|ψ>
    mul!(dstate, eq.imag_evo, state, -0.5, one(t))
    return
end

function Base.show(io::IO, mime::MIME"text/plain", eq::NoisySchrodingerEquation) #Pretty printing
    Base.show(io, mime, eq.equation)
    printstyled(io, "collapse operators: ", length(eq.num_cops); color=:yellow)
    return println(io)
end

"""
    struct NoisySchrodingerEquation
    NoisySchrodingerEquation(reg, tspan, hamiltonian, c_ops; kw...)
    NoisySchrodingerEquation(reg, tspan, hamiltonian, noise_model; kw...)

Define an ODE Problem representing a single noisy trajectory in an open system in the weak-coupling limit

# Arguments

- `reg`: evolution register and initial state of the problem
- `tspan`: a vector of times at which to save the simulation 
- `hamiltonian`: AbstractBlock such as that produced by `rydberg_h`` 
    representing evolution hamiltonian
"""
struct NoisySchrodingerProblem{Reg,EquationType<:ODEFunction,uType,tType,Algo,Kwargs} <:
       SciMLBase.AbstractODEProblem{uType,tType,true}
    reg::Reg
    f::EquationType
    c_ops::Vector
    state::uType
    u0::uType
    tspan::tType
    algo::Algo
    coherent_noise::Function
    confusion_mat #Matrix-like
    save_times
    kwargs::Kwargs
    p::Real
end

#internal constructor with all of the options
function NoisySchrodingerProblem(
    reg::AbstractRegister,
    save_times,
    expr,
    c_ops::Vector,
    coherent_noise::Function,
    confusion_mat; #TODO: property checking. Needs to support matrix-vector multiplication
    algo=DP8(), #Algorithm that was used for testing
    kw... #proceed with caution, some kwargs might break things
)
    nqudits(reg) == nqudits(expr) || throw(ArgumentError("number of qubits/sites does not match!"))

    issorted(save_times) || throw(ArgumentError("save_times must be sorted"))
    tspan = (0.0, last(save_times)) #Always starting at t=0 to avoid silly mistakes

    state = statevec(reg) #get initial statevector
    space = YaoSubspaceArrayReg.space(reg) #BloqadeNoisy does not currently support subspaces

    #type check c_ops
    try
        if length(c_ops) != 0
            [c^2 * state for c in c_ops] #making sure matrix-vector product works and matrix is square
        else
            c_ops = [0 * I]
        end
    catch
        throw(ArgumentError("Collapse operators are either not square or do not match register type"))
    end

    #type check confusion_mat
    try
        confusion_mat^2 * abs.(state) .^ 2 #check that matrix-vector product works and matrix is square
    catch
        throw(ArgumentError("Confusion matrix does not match register type"))
    end

    #type check coherent noise
    test_hamiltonian_sample = coherent_noise()
    typeof(test_hamiltonian_sample) == typeof(expr) || throw(ArgumentError("Coherent noise sample is not compatible with Hamiltonian"))
    nqudits(test_hamiltonian_sample) == nqudits(expr) || throw(ArgumentError("Coherent noise sample is not compatible with Hamiltonian"))

    T = real(eltype(state))
    T = isreal(expr) ? T : Complex{T}
    eq = NoisySchrodingerEquation(
        SchrodingerEquation(
            expr, Hamiltonian(T, expr, space)
        ),
        sum([l' * l for l in c_ops]),
        length(c_ops)
    )
    ode_f = ODEFunction(eq)

    tspan_type = promote_type(real(eltype(state)), eltype(tspan))
    tspan = tspan_type.(tspan) # promote tspan to T so Dual number works

    jump_cb = ContinuousCallback( #trigger quantum jumps and collapse the state
        collapse_condition,
        (integrator) -> collapse!(integrator, c_ops),
    )

    #remove save_start and sate_on options to allow saving at specified times
    ode_options = (
        save_everystep=false, #save only at designated save_times, improves performance
        save_start=true,
        dense=false,
        reltol=1e-10,
        abstol=1e-10,
        saveat=save_times,
        callback=jump_cb #quantum jumps
    )

    kw = pairs(merge(ode_options, kw))

    return NoisySchrodingerProblem(
        reg,
        ode_f,
        c_ops,
        state,
        copy(state),
        tspan,
        algo,
        coherent_noise,
        confusion_mat,
        save_times,
        kw,
        rand(), #random initial condition
    )
end

function NoisySchrodingerProblem( #Generic simulation with collapse operators
    reg::ArrayReg,
    save_times,
    h,
    c_ops::Array;
    kw...
)
    return NoisySchrodingerProblem(
        reg,
        save_times,
        h,
        c_ops,
        () -> h,
        Matrix(I, 2^nqubits(reg), 2^nqubits(reg));
        kw...
    )
end

function NoisySchrodingerProblem( #Simulation based on NoiseModel
    reg::ArrayReg,
    save_times,
    h::AbstractBlock,
    noise_model::ErrorModel;
    kw...
)
    return NoisySchrodingerProblem(
        reg,
        save_times,
        h,
        noise_model.collapse_ops(nqubits(reg)),
        noise_model.coherent_noise(h),
        noise_model.confusion_mat(nqubits(reg));
        kw...
    )
end

function collapse_condition(u, t, integrator) #condition for continuous callback
    norm(u)^2 - integrator.p
end

function collapse!(integrator, L_ops) #trigger a quantum jump
    dp = 0
    l = length(L_ops)
    probs = Vector{Float64}(undef, l)
    for i in 1:l
        dp += norm(L_ops[i] * integrator.u)^2 #normalization
        probs[i] = dp #cumulative distribution
    end
    r = rand()
    for i in 1:l
        if r <= probs[i] / dp #choose jump based on r
            copy!(integrator.u, L_ops[i] * integrator.u) #jump
            normalize!(integrator.u) #normalize
            break
        end
    end
    integrator.p = rand() #randomize for next jump
end

"""
    function randomize

Used to randomly modify the parameters of the Hamiltonian and randomly
change the initial condition, producing another NoisySchrodingerProblem object
representing a new trajectory. This is a replacement for SciML remake
"""
function randomize(prob::NoisySchrodingerProblem)
    h = prob.coherent_noise()
    p = rand()
    space = YaoSubspaceArrayReg.space(prob.reg)
    T = real(eltype(prob.state))
    T = isreal(h) ? T : Complex{T}
    eq = NoisySchrodingerEquation(
        SchrodingerEquation(
            h, Hamiltonian(T, h, space)
        ),
        sum([l' * l for l in prob.c_ops]),
        length(prob.c_ops)
    )
    ode_f = ODEFunction(eq)
    return NoisySchrodingerProblem(
        copy(prob.reg),
        ode_f,
        prob.c_ops,
        copy(prob.state),
        copy(prob.u0),
        prob.tspan,
        prob.algo,
        prob.coherent_noise,
        prob.confusion_mat,
        prob.save_times,
        prob.kwargs,
        p
    )
end

function Base.show(io::IO, mime::MIME"text/plain", prob::NoisySchrodingerProblem) #pretty printing
    indent = get(io, :indent, 0)
    tab(indent) = " "^indent

    println(io, tab(indent), "NoisySchrodingerProblem:")
    # state info
    println(io, tab(indent + 2), "register info:")
    print(io, tab(indent + 4), "type: ")
    printstyled(io, typeof(prob.reg); color=:green)
    println(io)

    print(io, tab(indent + 4), "storage size: ")
    printstyled(io, Base.format_bytes(storage_size(prob.reg)); color=:yellow)
    println(io)
    println(io)

    # tspan info
    println(io, tab(indent + 2), "time span (μs): ", prob.tspan)
    println(io)

    # equation info
    println(io, tab(indent + 2), "equation: ")
    show(IOContext(io, :indent => indent + 4), mime, prob.f.f)
    println(io, tab(indent + 4), "algorithm: ", repr(prob.algo))
end

"""
    function solve

This method should generally not be called directly unless you are implementing your own EnsembleProblem. This Implements
the DifferentialEquations interface for NoisySchrodingerProblem, which represents a single quantum trajectory.
"""
function DiffEqBase.solve(prob::NoisySchrodingerProblem, args...; sensealg=nothing, initial_state=nothing, kw...)
    if sensealg === nothing && haskey(prob.kwargs, :sensealg)
        sensealg = prob.kwargs[:sensealg]
    end
    # update initial state
    if initial_state !== nothing
        initial_state isa AbstractRegister ||
            throw(ArgumentError("initial_state must be a register, got $(typeof(initial_state))"))
        u0 = statevec(initial_state)
    else
        u0 = prob.u0
    end
    return DiffEqBase.solve_up(prob, sensealg, u0, nothing, args...; kw...)
end

DiffEqBase.get_concrete_problem(prob::NoisySchrodingerProblem, isadapt; kw...) = prob

"""
    function emulate

Emulate the evolution of a noisy system by averaging over an ensemble of ntraj noisy trajectories. The output is an array
where each array contains the output of output_func on the timeseries corresponding to a single trajectory.

# Arguments

- `prob`: NoisySchrodingerProblem to emulate
- `ntraj`: number of trajectories to use for simulation
- `output_func`: Vector{Complex}->Any - Transformation of the statevector array to save
- `ensemble_algo`: See EnsembleProblem documentation. Common choices are EnsembleSerial and EnsembleThreaded
"""
function emulate_noisy(
    prob::NoisySchrodingerProblem,
    ntraj::Int,
    output_func::Function;
    ensemble_algo=EnsembleSerial()
)
    ensemble_prob = EnsembleProblem(
        prob,
        prob_func=(prob, i, repeat) -> randomize(prob),
        output_func=(sol, i) -> (output_func([normalize(sol(t)) for t in prob.save_times]), false)
    )
    solve(ensemble_prob, prob.algo, ensemble_algo, trajectories=ntraj)
end

"""
    function emulate

Emulate an ensemble and return the ensemble average over bitstring distribution.

# Arguments
- `prob`: NoisySchrodingerProblem
- `ntraj`: Number of trajectories
- `report_error`: Returns 2σ estimated from the sample
- `ensemble_algo`: What type of parallelization is desired. See EnsembleSimulation documentation
"""
function emulate_noisy(
    prob::NoisySchrodingerProblem,
    ntraj::Int;
    report_error=false,
    ensemble_algo=EnsembleSerial()
)

    output_func = sol -> [abs.(u) .^ 2 for u in sol] #store the probability distribution
    sim = emulate_noisy(prob, ntraj, output_func; ensemble_algo=ensemble_algo) #call up
    return if report_error
        (
            amps=simulation_series_mean(sim),
            twosigma=simulation_series_err(sim)
        )
    else
        simulation_series_mean(sim)
    end
end

"""
    function emulate

Emulate an ensemble and return the ensemble average of a set of expectation values

# Arguments
- `prob`: NoisySchrodingerProblem
- `ntraj`: Number of trajectories
- `expectations`: Set of Hermitian operators representing observables. Matrix-like objects
- `shots`: Whether to resample the solution at a fixed number of shots. No resamping if `shots = 0``
- `report_error`: Returns 2σ estimated from the sample. Requires expectation values to be diagonal in computational basis
- `ensemble_algo`: What type of parallelization is desired. See EnsembleSimulation documentation
"""
function emulate_noisy(
    prob::NoisySchrodingerProblem,
    ntraj::Int,
    expectations::Array;
    readout_error=false,
    shots=0,
    report_error=false,
    ensemble_algo=EnsembleSerial()
)

    @assert all(ishermitian.(expectations)) #observables are hermitian!

    if !readout_error #Was not sure how to do this cleanly because multiple dispatch does not account for different kwargs...
        output_func = (sol) -> [[real(u' * (e * u)) for e in expectations] for u in sol]
        sim = emulate_noisy(prob, ntraj, output_func; ensemble_algo=ensemble_algo)

        return if report_error
            (
                expectations=[simulation_series_mean(sim; index = i) for (i, _) in enumerate(expectations)],
                twosigma=[simulation_series_err(sim; index = i) for (i, _) in enumerate(expectations)]
            )
        else
            [simulation_series_mean(sim; index=i) for (i, _) in enumerate(expectations)]
        end
    else
        @assert all([typeof(e) <: Diagonal for e in expectations])

        sim = emulate_noisy(prob, ntraj; report_error=true, ensemble_algo=ensemble_algo)
        amps = sim.amps
        amps_err = sim.twosigma

        if shots == 0 #no resampling
            res = [[_expectation_value_noisy(prob.confusion_mat, a, e, err) for (a, err) in zip(amps, amps_err)] for e in expectations]
            expec, perr = [[[a[i] for a in e] for e in res] for i in 1:2]
            return if report_error
                (
                    expectations=expec,
                    propagated_error=perr
                )
            else
                return expec
            end
        else
            res = [[_expectation_value_noisy(prob.confusion_mat, a, e, err, shots) for (a, err) in zip(amps, amps_err)] for e in expectations]
            mval, sample_err, err = [[[a[i] for a in e] for e in res] for i in 1:3]
            return if report_error
                (
                    expectations=mval,
                    shot_error=sample_err,
                    propagated_err=err
                )
            else
                return mval
            end
        end
    end
end

"""
    function expec_series_mean

Convenience method to access the expectation value over the ensemble
over the series of save times.

# Arguments
- sim: EnsembleSolution, result of calling `emulate`
- index: index of the expectation value in the array provided
"""
function simulation_series_mean(sim; index=false)
    ntraj = length(sim)
    times = length(sim[1])
    if index == false
        [mean([sim[i][t] for i in 1:ntraj]) for t in 1:times]
    else
        [mean([sim[i][t][index] for i in 1:ntraj]) for t in 1:times]
    end
end

"""
    function expec_series_err

Convenience method to estimate the sampling error in the ensemble solution

# Arguments
- `sim`: EnsembleSolution, result of calling `emulate`
- `index`: index of the expectation value in the array provided
- `factor`: How many σ from the mean to report
"""
function simulation_series_err(sim; index=false, factor=2)
    ntraj = length(sim)
    times = length(sim[1])
    if index == false
        [factor * std([sim[i][t] for i in 1:ntraj]) / sqrt(ntraj) for t in 1:times]
    else
        [factor * std([sim[i][t][index] for i in 1:ntraj]) / sqrt(ntraj) for t in 1:times]
    end
end

function _expectation_value_noisy(
    cmat,
    amps,
    op::Diagonal,
    errs::Vector
)
    expec = sum([a * real(n) for (a, n) in zip(cmat * amps, op.diag)])
    (
        expec,
        sqrt(sum([(err * real(n))^2 for (err, n) in zip(cmat * errs, op.diag)]))
    )
end

function _expectation_value_noisy(
    cmat,
    amps,
    op::Diagonal,
    errs::Vector,
    shots::Int
)
    w = Weights(cmat * amps) #create weights representing measurement probabilities
    samples = sample(1:length(amps), w, shots)
    S = [real(op.diag[s]) for s in samples]
    (
        mean(S),
        2 * std(S) / sqrt(shots), #Using 2σ because 95% of the data is within this error
        sqrt(sum([(err * real(n))^2 for (err, n) in zip(cmat * errs, op.diag)]))
    )
end