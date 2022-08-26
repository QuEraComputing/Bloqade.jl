using OnlineStats


abstract type AbstractRunStats; end

struct NoStats <: AbstractRunStats; end
@inline OnlineStats.fit!(R::NoStats, ::Symbol, val) = R


struct RunStats{T <: Real} <: AbstractRunStats
    diag_update_fails::Variance{T}
    cluster_update_accept::Variance{T}
    cluster_count::Variance{T}
    cluster_sizes::Variance{T}
    accepted_cluster_count::Variance{T}
    accepted_cluster_sizes::Variance{T}
    rejected_cluster_count::Variance{T}
    rejected_cluster_sizes::Variance{T}

    RunStats{T}() where T = new{T}(Variance(T), Variance(T), Variance(T), Variance(T),
                                   Variance(T), Variance(T), Variance(T), Variance(T))
end
RunStats() = RunStats{Float64}()
@inline OnlineStats.fit!(R::RunStats{T}, field::Symbol, val) where T = (fit!(getproperty(R, field), T(val)); R)


struct RunStatsHistogram{T <: Real, R <: StepRangeLen} <: AbstractRunStats
    diag_update_fails::Hist{T, R}
    cluster_update_accept::Hist{T, R}
    cluster_count::Vector{Int}
    cluster_sizes::Vector{Int}
    accepted_cluster_count::Vector{Int}
    accepted_cluster_sizes::Vector{Int}
    rejected_cluster_count::Vector{Int}
    rejected_cluster_sizes::Vector{Int}

    RunStatsHistogram{T}(size::Int) where T =
        new{T, typeof(0:(1/size):1)}(Hist(0:(1/size):1, T), Hist(0:(1/size):1, T),
                                     zeros(Int, size), zeros(Int, size),
                                     zeros(Int, size), zeros(Int, size),
                                     zeros(Int, size), zeros(Int, size))
end
RunStatsHistogram(size::Int) = RunStatsHistogram{Float64}(size::Int)
function OnlineStats.fit!(R::RunStatsHistogram{T}, field::Symbol, val) where T
    if field == :diag_update_fails || field == :cluster_update_accept
        fit!(getproperty(R, field), T(val))
    else
        V = getfield(R, field)
        old_len = length(V)
        if val > old_len
            Base.resize!(V, val)
            V[(old_len+1):end] .= 0
        end
        V[val] += 1
    end
    return R
end
