function draw(wf::Waveform)
    clocks = sample_clock(wf)
    fig = plt.plot(
        clocks, EaRydWaveforms._rm_err.(sample_values(wf, clocks));
    )
    plt.xlabel("clock (μs)")
    plt.ylabel("value (rad/µs)")
    return fig
end
