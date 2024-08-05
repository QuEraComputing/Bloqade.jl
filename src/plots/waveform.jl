function plot!(ax, wf::Waveform; kwargs...)
    clocks = sample_clock(wf)
    lines!(ax, clocks, BloqadeWaveforms._rm_err.(sample_values(wf, clocks) ./ (2π)); kwargs...)
    return ax
end

function plot(wf::Waveform; kwargs...)
    fig = CairoMakie.Figure()
    ax = Axis(fig[1, 1]; xlabel="time (μs)", ylabel = "value / 2π (MHz)")
    plot!(ax, wf; kwargs...)
    return fig
end