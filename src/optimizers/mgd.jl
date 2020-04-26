export mgd

# prirating!
function LsqFit.check_data_health(xdata::Vector{<:AbstractVector}, ydata)
    true
end

"""
    mgd(f, x0; γ, δ, k, α, A, ξ, ϵ, n, p0)

Model gradient descent. `p0` is the initial parameters of the quadratic model, e.g. if your quadratic models has two variables, `p0` should be a vector of length `6`.
For other parameters, see the appendix of: https://arxiv.org/pdf/2004.04197.pdf
"""
function mgd(f, x0::XT; γ, δ, k, α, A, ξ, ϵ, n, p0) where XT
    m = 0
    x = x0
    params = p0
    nx = length(x)
    Lx = XT[]
    Ly = Float64[]
    neval = 0
    while neval + k <= n
        push!(Lx, copy(x))
        push!(Ly, f(x))
        println("Step $m, loss = $(Ly[end])")
        δi = δ/(m+1)^ξ
        for i in 1:k
            x2 = rand_disk(x, δi)
            push!(Lx, x2)
            push!(Ly, f(x2))
            neval += 1
        end
        Lx2 = XT[]
        Ly2 = Float64[]
        for (x2, y2) in zip(Lx, Ly)
            if norm(x2 - x) <= δi
                push!(Lx2, x2)
                push!(Ly2, y2)
            end
        end
        params = _fit(multivariate_quadratic, Lx2, Ly2, params)
        g = NiLang.AD.gradient(Val(1), multivariate_quadratic, (0.0, x, params))[2]
        γ2 = γ/(m+1+A)^α
        if γ2*norm(g) < ϵ
            return x
        end
        x .-= γ2 .* g
        m += 1
    end
    return x
end

function rand_disk(x, δi)
    while true
        np = rand(length(x)) .* 2δi
        if norm(np) < δi
            return x .+ np
        end
    end
end

function _fit(fit_func, xs, ys, p0)
    function model(xs, p)
        map(x->fit_func(0.0, x, p)[1], xs)
    end
    fit = curve_fit(model, xs, ys, p0)
    println("fitting error = $(fit.resid[end])")
    fit.param
end

@i function multivariate_quadratic(res::T, x::AbstractVector, p) where T
    k ← 1
    res += identity(p[k])
    @invcheckoff @inbounds for j = 1:length(x)
        k += identity(1)
        res += x[j] * p[k]
    end
    @invcheckoff @inbounds for j = 1:length(x)
        anc ← zero(T)
        for i=1:j-1
            k += identity(1)
            anc += x[i] * x[j]
            res += anc * p[k]
            anc -= x[i] * x[j]
        end
        k += identity(1)
        anc += x[j] ^ 2
        res +=  anc * p[k]
        anc -= x[j] ^ 2
    end
    k → length(x)*(length(x)+1) ÷ 2 + length(x) + 1
end
