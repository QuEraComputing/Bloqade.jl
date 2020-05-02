"""
Spall's implementation of SPSA algorithm.

    Spall, J. C. (1998).
    Implementation of the simultaneous perturbation algorithm for stochastic optimization.
    IEEE Transactions on Aerospace and Electronic Systems.
    https://doi.org/10.1109/7.705889
"""

export spsa

"""
second-order simultaneous perturbation stochastic approximation.

Args:
    f: loss function.
    x0: initial variables.
    bounds: lower bound and higher bound for variables, None for no bounds.
    ac: initial learning rate and initial perturbation stength.
    A: statbility constant.
    alpha: decay rate for learning speed.
    gamma: decay rate for perturbation strength.
    maxiter: maximum number of iteration.

Note:
    The choice of parameters,
    * `alpha` and `gamma` have thoretically valid default values 0.602 and 0.101.
    * in hight noise setting, pick smaller `a` and larger `c`.
"""
function spsa(f, x0::AbstractVector{T};  bounds=nothing, ac=(0.2, 0.5), alpha=0.602,
          A=nothing, gamma=0.101, maxiter=5000, secondorder=false) where T
    if A === nothing
        A = 0.1 * maxiter
    end
    a, c = ac

    if secondorder
        Hk = zeros(T, length(x0), length(x0))
    end
    for k in 1:maxiter
        ak = a / (k + A)^alpha
        ck = c / k^gamma
        if secondorder
            g, delta, pos, neg = _get_g(f, x0, ck)
            Hk_ = _get_h(f, x0, ck, delta, pos, neg)
            Hk .= ((k-1) / k) .* Hk .+ (1 / k) .* Hk_
            x0 .-= ak .* (regularted_inv(Hk; delta=1e-5) * g)
        else
            g, fpos, fneg = _get_g(f, x0, ck)
            x0 .-= ak .* g
        end
        @show f(x0)
        if bounds !== nothing
            clip!(x0, bounds...)
        end
    end
    #return OptimizeResult(x=x0, fun=f(x0), success=true)
    return x0
end

function clip!(x, a, b)
    @inbounds for i = eachindex(x)
        if x[i] < a
            x[i] = a
        elseif x[i] > b
            x[i] = b
        end
    end
end

"""
suggest hyper parameters for SPSA.

Args:
    fun (func): loss function.
    x0 (ndarray): initial variables.
    initial_step (float): the initial step size.
    noise_strength (float): the noise level (standard error of fun).

Return:
    (a, c): hyper parameters.
"""
function _hp_suggest(fun, x0, initial_step, noise_strength, alpha, A)
    num_eval = 10  # evaluate gradient 10 times to get the mean step.
    c = noise_strength * 300 + 1e-4
    g0 = sum(x->_get_g(fun, x0, c, return_hessian=false), 1:num_eval)/num_eval
    # use g0*ak = initial_step to set parameters
    a = initial_step / g0 * (1 + A)^alpha
    return a, c
end


"""
compute the gradient
"""
function _get_g(fun, x0, ck)
    p = length(x0)
    delta = rand([-1,1], p) .* ck
    xpos = x0 .+ delta
    xneg = x0 .- delta
    fpos, fneg = fun(xpos), fun(xneg)
    (fpos - fneg) ./ (2 .* delta), delta, (xpos, fpos), (xneg, fneg)
end

function _get_h(fun, x0, ck, delta, pos, neg)
    xpos, fpos = pos
    xneg, fneg = neg
    delta1 = rand([-1,1], length(x0)) .* ck
    fneg_ = fun(xneg .+ delta1)
    fpos_ = fun(xpos .+ delta1)
    g1n = (fneg_ - fneg) ./ delta1
    g1p = (fpos_ - fpos) ./ delta1
    hessian = (g1p .- g1n) ./ 4. ./ delta'
    hessian .+= hessian'
    return hessian
end

function regularted_inv(A; delta)
    inv(sqrt(A * A) + delta * I)
end
