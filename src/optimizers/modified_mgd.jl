"""
    modified_mgd(f, x0; δ, k, ϵ, n, p0, lr=0.5)

Modified model gradient descent.
See `mgd` for input parameter informations.
"""
function modified_mgd(f, x0::XT; δ, k, ϵ, n, p0, lr=0.5) where XT
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
        for i in 1:k
            x2 = rand_disk(x, δ)
            push!(Lx, x2)
            push!(Ly, f(x2))
            neval += 1
        end
        Lx2 = XT[]
        Ly2 = Float64[]
        for (x2, y2) in zip(Lx, Ly)
            if norm(x2 - x) <= δ + 1e-5
                push!(Lx2, x2)
                push!(Ly2, y2)
            end
        end
        params = _fit(multivariate_quadratic, Lx2, Ly2, params)
        x2, ismin = _optimal(multivariate_quadratic, params)
        if norm(x .- x2) < ϵ
            return x2
        else
            δ = norm(x .- x2)
            x = ismin ? x + (x2-x)*lr : x - (x2-x)*lr
        end
        m += 1
    end
    return x
end

function _optimal(::typeof(multivariate_quadratic), params::AbstractVector{T}) where T
    nx = (isqrt(8*length(params) + 1) - 3) ÷ 2
    b = -params[2:1+nx]
    A = zeros(T, nx, nx)
    k = 1+nx
    for i=1:nx
        for j=1:i-1
            k += 1
            A[i,j] = params[k]
            A[j,i] = params[k]
        end
        k += 1
        A[i,i] = 2*params[k]
    end
    if isposdef(A)
        ismin = true
    else
        ismin = false
    end
    gmres(A, b), ismin
end
