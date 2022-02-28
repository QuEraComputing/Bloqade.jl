function draw(wf::Waveform; kw...)
    fig = Figure(resolution=(1400, 600))
    draw!(fig[1, 1], wf; kw...)
    return fig
end

function draw!(plt, wf::Waveform; title::String="")
    clocks = sample_clock(wf)
    ax = Axis(
        plt;
        # xticks = (1:length(indices), string.(indices.-1;base=2, pad=nqubits(r))),
        title,
        xlabel="clock (μs)",
        ylabel="value (rad/µs)",
    )
    lines!(ax, clocks, EaRydWaveforms._rm_err.(sample_values(wf, clocks)))
    return plt
end
