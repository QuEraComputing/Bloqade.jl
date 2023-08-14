# this function is copied from https://github.com/SciML/OrdinaryDiffEq.jl/blob/HEAD/src/solve.jl#L9-L470
# should be replaced with a cleaner implementation when __init gets improved in the future
#
# this function has the following modification:
# 1. we need to let integrator modify our own preallocated state member instead of u0
# 2. we need to get rid of the interpolation
# 3. kwargs alias_du0 and alias_u0 are removed, since we don't need them anymore

# TODOs
# we should have the default options directly be our preferred (save_everystep=false, etc.)
# but this should just be inside the problem type instead of the init function.
using Pkg
get_ver = pkg_name -> filter(x-> x.second.name == pkg_name, Pkg.dependencies()) |> x -> first(x)[2].version   


vname = :(destats)
mname = :(DEStats)
begin
    if get_ver("DiffEqBase") >= v"6.122.0"
        vname = :(stats)
        mname = :(Stats)
    end
end


function DiffEqBase.__init(
    prob::SchrodingerProblem,
    alg::OrdinaryDiffEqAlgorithm,
    timeseries_init = (),
    ts_init = (),
    ks_init = (),
    recompile::Type{Val{recompile_flag}} = Val{true};
    saveat = (),
    tstops = (),
    d_discontinuities = (),
    save_idxs = nothing,
    save_everystep = false,
    save_on = false,
    save_start = false,
    save_end = nothing,
    callback = nothing,
    dense = save_everystep && !(typeof(alg) <: Union{DAEAlgorithm,FunctionMap}) && isempty(saveat),
    calck = (callback !== nothing && callback !== CallbackSet()) || (dense) || !isempty(saveat), # and no dense output
    dt = alg isa FunctionMap && isempty(tstops) ? eltype(prob.tspan)(1) : eltype(prob.tspan)(0),
    dtmin = nothing,
    dtmax = eltype(prob.tspan)((prob.tspan[end] - prob.tspan[1])),
    force_dtmin = false,
    adaptive = isadaptive(alg),
    gamma = gamma_default(alg),
    abstol = nothing,
    reltol = nothing,
    qmin = qmin_default(alg),
    qmax = qmax_default(alg),
    qsteady_min = qsteady_min_default(alg),
    qsteady_max = qsteady_max_default(alg),
    beta1 = nothing,
    beta2 = nothing,
    qoldinit = isadaptive(alg) ? 1 // 10^4 : 0,
    controller = nothing,
    fullnormalize = true,
    failfactor = 2,
    maxiters = adaptive ? 1000000 : typemax(Int),
    internalnorm = ODE_DEFAULT_NORM,
    internalopnorm = LinearAlgebra.opnorm,
    isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
    unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
    verbose = true,
    timeseries_errors = true,
    dense_errors = false,
    advance_to_tstop = false,
    stop_at_next_tstop = false,
    initialize_save = true,
    progress = false,
    progress_steps = 1000,
    progress_name = "ODE",
    progress_message = ODE_DEFAULT_PROG_MESSAGE,
    userdata = nothing,
    allow_extrapolation = alg_extrapolates(alg),
    initialize_integrator = true,
    initializealg = DefaultInit(),
    kwargs...,
) where {recompile_flag}
    if typeof(prob.f) <: DynamicalODEFunction && typeof(prob.f.mass_matrix) <: Tuple
        if any(mm != I for mm in prob.f.mass_matrix)
            error("This solver is not able to use mass matrices.")
        end
    elseif !(typeof(prob) <: DiscreteProblem) &&
           !(typeof(prob) <: DiffEqBase.AbstractDAEProblem) &&
           !is_mass_matrix_alg(alg) &&
           prob.f.mass_matrix != I
        error("This solver is not able to use mass matrices.")
    end

    if !isempty(saveat) && dense
        @warn(
            "Dense output is incompatible with saveat. Please use the SavingCallback from the Callback Library to mix the two behaviors."
        )
    end

    progress && @logmsg(LogLevel(-1), progress_name, _id = _id = :OrdinaryDiffEq, progress = 0)

    tType = eltype(prob.tspan)
    tspan = prob.tspan
    tdir = sign(tspan[end] - tspan[1])

    t = tspan[1]

    if (
        (
            (
                !(typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm) &&
                !(typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm) &&
                !(typeof(alg) <: DAEAlgorithm)
            ) || !adaptive
        ) &&
        dt == tType(0) &&
        isempty(tstops)
    ) && !(typeof(alg) <: Union{FunctionMap,LinearExponential})
        error("Fixed timestep methods require a choice of dt or choosing the tstops")
    end

    isdae =
        alg isa DAEAlgorithm || (
            !(typeof(prob) <: DiscreteProblem) &&
            prob.f.mass_matrix != I &&
            !(typeof(prob.f.mass_matrix) <: Tuple) &&
            ArrayInterface.issingular(prob.f.mass_matrix)
        )
    if alg isa CompositeAlgorithm && alg.choice_function isa AutoSwitch
        auto = alg.choice_function
        _alg = CompositeAlgorithm(
            alg.algs,
            AutoSwitchCache(
                0,
                0,
                auto.nonstiffalg,
                auto.stiffalg,
                auto.stiffalgfirst,
                auto.maxstiffstep,
                auto.maxnonstiffstep,
                auto.nonstifftol,
                auto.stifftol,
                auto.dtfac,
                auto.stiffalgfirst,
                auto.switch_max,
            ),
        )
    else
        _alg = alg
    end
    f = prob.f
    p = prob.p

    # Get the control variables

    u = copyto!(prob.state, prob.u0)
    du = nothing
    duprev = nothing

    uType = typeof(u)
    uBottomEltype = recursive_bottom_eltype(u)
    uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

    uEltypeNoUnits = recursive_unitless_eltype(u)
    tTypeNoUnits = typeof(one(tType))

    if typeof(_alg) <: FunctionMap
        abstol_internal = false
    elseif abstol === nothing
        if uBottomEltypeNoUnits == uBottomEltype
            abstol_internal = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1 // 10^6))
        else
            abstol_internal = real.(oneunit.(u) .* 1 // 10^6)
        end
    else
        abstol_internal = real.(abstol)
    end

    if typeof(_alg) <: FunctionMap
        reltol_internal = false
    elseif reltol === nothing
        if uBottomEltypeNoUnits == uBottomEltype
            reltol_internal = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1 // 10^3))
        else
            reltol_internal = real.(oneunit.(u) .* 1 // 10^3)
        end
    else
        reltol_internal = real.(reltol)
    end

    dtmax > zero(dtmax) && tdir < 0 && (dtmax *= tdir) # Allow positive dtmax, but auto-convert
    # dtmin is all abs => does not care about sign already.

    if !isdae &&
       isinplace(prob) &&
       typeof(u) <: AbstractArray &&
       eltype(u) <: Number &&
       uBottomEltypeNoUnits == uBottomEltype &&
       tType == tTypeNoUnits # Could this be more efficient for other arrays?
        rate_prototype = recursivecopy(u)
    elseif prob isa DAEProblem
        rate_prototype = prob.du0
    else
        if (uBottomEltypeNoUnits == uBottomEltype && tType == tTypeNoUnits) || eltype(u) <: Enum
            rate_prototype = u
        else # has units!
            rate_prototype = u / oneunit(tType)
        end
    end
    rateType = typeof(rate_prototype) ## Can be different if united

    if isdae
        if uBottomEltype == uBottomEltypeNoUnits
            res_prototype = u
        else
            res_prototype = one(u)
        end
        resType = typeof(res_prototype)
    end

    tstops_internal = initialize_tstops(tType, tstops, d_discontinuities, tspan)
    saveat_internal = initialize_saveat(tType, saveat, tspan)
    d_discontinuities_internal = initialize_d_discontinuities(tType, d_discontinuities, tspan)

    callbacks_internal = CallbackSet(callback)

    max_len_cb = DiffEqBase.max_vector_callback_length_int(callbacks_internal)
    if max_len_cb !== nothing
        uBottomEltypeReal = real(uBottomEltype)
        if isinplace(prob)
            callback_cache = DiffEqBase.CallbackCache(u, max_len_cb, uBottomEltypeReal, uBottomEltypeReal)
        else
            callback_cache = DiffEqBase.CallbackCache(max_len_cb, uBottomEltypeReal, uBottomEltypeReal)
        end
    else
        callback_cache = nothing
    end

    ### Algorithm-specific defaults ###
    if save_idxs === nothing
        ksEltype = Vector{rateType}
    else
        ks_prototype = rate_prototype[save_idxs]
        ksEltype = Vector{typeof(ks_prototype)}
    end

    # Have to convert incase passed in wrong.
    if save_idxs === nothing
        timeseries = timeseries_init === () ? uType[] : convert(Vector{uType}, timeseries_init)
    else
        u_initial = u[save_idxs]
        timeseries = timeseries_init === () ? typeof(u_initial)[] : convert(Vector{uType}, timeseries_init)
    end

    ts = ts_init === () ? tType[] : convert(Vector{tType}, ts_init)
    ks = ks_init === () ? ksEltype[] : convert(Vector{ksEltype}, ks_init)
    alg_choice = typeof(_alg) <: CompositeAlgorithm ? Int[] : ()

    if !adaptive && save_everystep && tspan[2] - tspan[1] != Inf
        if dt == 0
            steps = length(tstops)
        else
            dtmin === nothing && (dtmin = DiffEqBase.prob2dtmin(prob; use_end_time = true))
            abs(dt) < dtmin && throw(ArgumentError("Supplied dt is smaller than dtmin"))
            steps = ceil(Int, internalnorm((tspan[2] - tspan[1]) / dt, tspan[1]))
        end
        sizehint!(timeseries, steps + 1)
        sizehint!(ts, steps + 1)
        sizehint!(ks, steps + 1)
    elseif save_everystep
        sizehint!(timeseries, 50)
        sizehint!(ts, 50)
        sizehint!(ks, 50)
    elseif !isempty(saveat_internal)
        sizehint!(timeseries, length(saveat_internal) + 1)
        sizehint!(ts, length(saveat_internal) + 1)
        sizehint!(ks, length(saveat_internal) + 1)
    else
        sizehint!(timeseries, 2)
        sizehint!(ts, 2)
        sizehint!(ks, 2)
    end

    QT, EEstT = if tTypeNoUnits <: Integer
        typeof(qmin), typeof(qmin)
    elseif prob isa DiscreteProblem
        # The QT fields are not used for DiscreteProblems
        constvalue(tTypeNoUnits), constvalue(tTypeNoUnits)
    else
        typeof(DiffEqBase.value(internalnorm(u, t))), typeof(internalnorm(u, t))
    end

    k = rateType[]

    if uses_uprev(_alg, adaptive) || calck
        uprev = recursivecopy(u)
    else
        # Some algorithms do not use `uprev` explicitly. In that case, we can save
        # some memory by aliasing `uprev = u`, e.g. for "2N" low storage methods.
        uprev = u
    end
    if allow_extrapolation
        uprev2 = recursivecopy(u)
    else
        uprev2 = uprev
    end

    if prob isa DAEProblem
        cache = alg_cache(
            _alg,
            du,
            u,
            res_prototype,
            rate_prototype,
            uEltypeNoUnits,
            uBottomEltypeNoUnits,
            tTypeNoUnits,
            uprev,
            uprev2,
            f,
            t,
            dt,
            reltol_internal,
            p,
            calck,
            Val(isinplace(prob)),
        )
    else
        cache = alg_cache(
            _alg,
            u,
            rate_prototype,
            uEltypeNoUnits,
            uBottomEltypeNoUnits,
            tTypeNoUnits,
            uprev,
            uprev2,
            f,
            t,
            dt,
            reltol_internal,
            p,
            calck,
            Val(isinplace(prob)),
        )
    end

    # Setting up the step size controller
    if (beta1 !== nothing || beta2 !== nothing) && controller !== nothing
        throw(
            ArgumentError(
                "Setting both the legacy PID parameters `beta1, beta2 = $((beta1, beta2))` and the `controller = $controller` is not allowed.",
            ),
        )
    end

    if (beta1 !== nothing || beta2 !== nothing)
        message = "Providing the legacy PID parameters `beta1, beta2` is deprecated. Use the keyword argument `controller` instead."
        Base.depwarn(message, :init)
        Base.depwarn(message, :solve)
    end

    if controller === nothing
        controller = default_controller(
            _alg,
            cache,
            qoldinit,
            beta1 === nothing ? nothing : beta1,
            beta2 === nothing ? nothing : beta2,
        )
    end

    dtmin === nothing && (dtmin = DiffEqBase.prob2dtmin(prob; use_end_time = false))

    save_end_user = save_end
    save_end =
        save_end === nothing ? save_everystep || isempty(saveat) || saveat isa Number || prob.tspan[2] in saveat :
        save_end

    opts = DEOptions{
        typeof(abstol_internal),
        typeof(reltol_internal),
        QT,
        tType,
        typeof(controller),
        typeof(internalnorm),
        typeof(internalopnorm),
        typeof(save_end_user),
        typeof(callbacks_internal),
        typeof(isoutofdomain),
        typeof(progress_message),
        typeof(unstable_check),
        typeof(tstops_internal),
        typeof(d_discontinuities_internal),
        typeof(userdata),
        typeof(save_idxs),
        typeof(maxiters),
        typeof(tstops),
        typeof(saveat),
        typeof(d_discontinuities),
    }(
        maxiters,
        save_everystep,
        adaptive,
        abstol_internal,
        reltol_internal,
        QT(gamma),
        QT(qmax),
        QT(qmin),
        QT(qsteady_max),
        QT(qsteady_min),
        QT(qoldinit),
        QT(failfactor),
        tType(dtmax),
        tType(dtmin),
        controller,
        internalnorm,
        internalopnorm,
        save_idxs,
        tstops_internal,
        saveat_internal,
        d_discontinuities_internal,
        tstops,
        saveat,
        d_discontinuities,
        userdata,
        progress,
        progress_steps,
        progress_name,
        progress_message,
        timeseries_errors,
        dense_errors,
        dense,
        save_on,
        save_start,
        save_end,
        save_end_user,
        callbacks_internal,
        isoutofdomain,
        unstable_check,
        verbose,
        calck,
        force_dtmin,
        advance_to_tstop,
        stop_at_next_tstop,
    )

    eval( :($vname = DiffEqBase.$mname(0)) )
    #stats = DiffEqBase.Stats(0)

    if typeof(_alg) <: OrdinaryDiffEqCompositeAlgorithm
        id = CompositeInterpolationData(f, timeseries, ts, ks, alg_choice, dense, cache)
        eval(
            quote
                sol = DiffEqBase.build_solution(
                    $prob,
                    $_alg,
                    $ts,
                    $timeseries,
                    dense = $dense,
                    k = $ks,
                    interp = $id,
                    alg_choice = $alg_choice,
                    calculate_error = false,
                    $vname = $vname,
                    #stats = stats,
                )
            end
        )
    else
        id = InterpolationData(f, timeseries, ts, ks, dense, cache)
        eval(
            quote
                sol = DiffEqBase.build_solution(
                    $prob,
                    $_alg,
                    $ts,
                    $timeseries,
                    dense = $dense,
                    k = $ks,
                    interp = $id,
                    calculate_error = false,
                    $vname = $vname,
                    #stats = stats,
                )
            end
        )
    end

    if recompile_flag == true
        FType = typeof(f)
        SolType = typeof(sol)
        cacheType = typeof(cache)
    else
        FType = Function
        if _alg isa OrdinaryDiffEqAlgorithm
            SolType = DiffEqBase.AbstractODESolution
            cacheType = OrdinaryDiffEqCache
        else
            SolType = DiffEqBase.AbstractDAESolution
            cacheType = DAECache
        end
    end

    # rate/state = (state/time)/state = 1/t units, internalnorm drops units
    # we don't want to differentiate through eigenvalue estimation
    eigen_est = inv(one(tType))
    tprev = t
    dtcache = tType(dt)
    dtpropose = tType(dt)
    iter = 0
    kshortsize = 0
    reeval_fsal = false
    u_modified = false
    EEst = EEstT(1)
    just_hit_tstop = false
    isout = false
    accept_step = false
    force_stepfail = false
    last_stepfail = false
    do_error_check = true
    event_last_time = 0
    vector_event_last_time = 1
    last_event_error = typeof(_alg) <: FunctionMap ? false : zero(uBottomEltypeNoUnits)
    dtchangeable = isdtchangeable(_alg)
    q11 = QT(1)
    success_iter = 0
    erracc = QT(1)
    dtacc = tType(1)
    reinitiailize = true
    saveiter = 0 # Starts at 0 so first save is at 1
    saveiter_dense = 0

    integrator = ODEIntegrator{
        typeof(_alg),
        isinplace(prob),
        uType,
        typeof(du),
        tType,
        typeof(p),
        typeof(eigen_est),
        typeof(EEst),
        QT,
        typeof(tdir),
        typeof(k),
        SolType,
        FType,
        cacheType,
        typeof(opts),
        fsal_typeof(_alg, rate_prototype),
        typeof(last_event_error),
        typeof(callback_cache),
        typeof(initializealg),
    }(
        sol,
        u,
        du,
        k,
        t,
        tType(dt),
        f,
        p,
        uprev,
        uprev2,
        duprev,
        tprev,
        _alg,
        dtcache,
        dtchangeable,
        dtpropose,
        tdir,
        eigen_est,
        EEst,
        QT(qoldinit),
        q11,
        erracc,
        dtacc,
        success_iter,
        iter,
        saveiter,
        saveiter_dense,
        cache,
        callback_cache,
        kshortsize,
        force_stepfail,
        last_stepfail,
        just_hit_tstop,
        do_error_check,
        event_last_time,
        vector_event_last_time,
        last_event_error,
        accept_step,
        isout,
        reeval_fsal,
        u_modified,
        reinitiailize,
        isdae,
        opts,
        eval(vname),
        #stats,
        initializealg,
    )

    if initialize_integrator
        if isdae
            DiffEqBase.initialize_dae!(integrator)
        end

        if save_start
            integrator.saveiter += 1 # Starts at 1 so first save is at 2
            integrator.saveiter_dense += 1
            copyat_or_push!(ts, 1, t)
            if save_idxs === nothing
                copyat_or_push!(timeseries, 1, integrator.u)
                copyat_or_push!(ks, 1, [rate_prototype])
            else
                copyat_or_push!(timeseries, 1, u_initial, Val{false})
                copyat_or_push!(ks, 1, [ks_prototype])
            end
        else
            saveiter = 0 # Starts at 0 so first save is at 1
            saveiter_dense = 0
        end

        initialize_callbacks!(integrator, initialize_save)
        initialize!(integrator, integrator.cache)

        if typeof(_alg) <: OrdinaryDiffEqCompositeAlgorithm && save_start
            # Loop to get all of the extra possible saves in callback initialization
            for i in 1:integrator.saveiter
                copyat_or_push!(alg_choice, i, integrator.cache.current)
            end
        end
    end

    handle_dt!(integrator)

    return integrator
end
