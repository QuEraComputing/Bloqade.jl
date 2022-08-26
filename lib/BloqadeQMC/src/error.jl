using Measurements

function mean_and_stderr(f::Function, x::Vector)
    y = f.(x)
    μ = mean(y)
    σ = std(y; mean=μ) / sqrt(length(x))

    return μ ± σ
end
mean_and_stderr(x::Vector) = mean_and_stderr(identity, x)


# jackknife estimates the error of a function of a mean (or several means), f(<x>)
function jackknife(f::Function, xs::Vector...)
    sum_xs = [sum(xs[i]) for i in 1:length(xs)]
    N = length(xs[1])

    f_J = zeros(N)
    @inbounds for i in 1:N
        xs_J = [(sum_xs[j] - xs[j][i]) / (N - 1) for j in eachindex(xs)]
        f_J[i] = f(xs_J...)
    end

    μ = mean(f_J)
    σ2 = (N - 1) * (mean(abs2, f_J) - μ^2)
    σ = σ2 < 0 ? 0.0 : sqrt(σ2)

    f_ = f((sum_xs / N)...)
    # bias = ((N - 1) * μ) - ((N - 1) * f_)

    μ′ = N*f_ - (N-1)*μ

    return μ′ ± σ
end
jackknife(xs::Vector...) = jackknife(identity, xs...)


function bootstrap(f::Function, xs::Vector...; nboot::Int=500)
    N = length(xs[1])

    f_B = zeros(nboot)

    for b in 1:nboot
        idx = rand(1:N, N)
        xs_b = [mean(@views x[idx]) for x in xs]
        f_B[b] = f(xs_b...)
    end

    μ = mean(f_B)
    σ2 = (N / (N - 1)) * (mean(abs2, f_B) - μ^2)
    σ = σ2 < 0 ? 0.0 : sqrt(σ2)

    f_ = f([mean(x) for x in xs]...)

    μ′ = 2*f_ - μ

    return μ′ ± σ
end
bootstrap(xs::Vector...; nboot::Int=500) = bootstrap(identity, xs...; nboot=nboot)
