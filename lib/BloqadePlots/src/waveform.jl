function draw!(ax, wf::Waveform)
    clocks = sample_clock(wf)
    fig = ax.plot(
        clocks, BloqadeWaveforms._rm_err.(sample_values(wf, clocks));
    )
    ax.set_xlabel("time (μs)")
    ax.set_ylabel("value (rad/µs)")
    return ax
end

function draw(wf::Waveform)
    fig, ax = plt.subplots()
    draw!(ax, wf)
    return fig
end
