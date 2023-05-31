function plot!(ax, wf::Waveform)
    clocks = sample_clock(wf)
    fig = ax.plot(clocks, BloqadeWaveforms._rm_err.(sample_values(wf, clocks) ./ (2π));)
    ax.set_xlabel("time (μs)")
    ax.set_ylabel("value / 2π (MHz)")
    return ax
end

function plot(wf::Waveform)
    fig, ax = plt.subplots()
    plot!(ax, wf)
    return fig
end
