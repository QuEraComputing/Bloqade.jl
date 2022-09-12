
using BloqadeWaveforms
using Roots: find_zero

# TODO: move this to BloqadeWaveforms
# function inserts ramps at the beginning and end to set correct initial and final values
function pin_waveform_edges(wf::Waveform,
    max_slope::Real,
    begin_value::Real,
    end_value::Real)

    duration = wf.duration



    t_begin = if !isapprox(wf(0.0), begin_value;atol=eps(),rtol=√eps())
        ramp_up =  (sign(wf(0.0)-begin_value)*max_slope)
        lin_ramp_begin = Waveform(t -> ramp_up .* t .+ begin_value, duration)
        find_zero(wf-lin_ramp_begin,(0.0,duration))
    else
        0.0
    end

    t_end = if !isapprox(wf(duration), end_value;atol=eps(),rtol=√eps())
        ramp_down = (sign(end_value-wf(duration))*max_slope)
        lin_ramp_end = Waveform(t -> ramp_down .* (t.-duration) .+ end_value, duration)
        find_zero(wf-lin_ramp_end,(0.0,duration))
    else
        duration
    end
        
    if t_begin > 0 && t_end < duration
        mid_wf = Waveform(t->wf.f(t.+t_begin),duration - t_begin - t_end)

        start_wf = linear_ramp(;
            duration=t_begin,
            start_value=begin_value,
            stop_value=wf(t_begin)
        )

        end_wf = linear_ramp(;
            duration=duration-t_end,
            start_value=wf(t_end),
            stop_value=end_value
        )

        return append(start_wf,mid_wf,end_wf)
    elseif t_begin > 0 
        end_wf = Waveform(t->wf.f(t.+t_begin),duration - t_begin)

        start_wf = linear_ramp(;duration=t_begin,
            start_value=begin_value,
            stop_value=wf(t_begin)
        )

        return append(start_wf,end_wf)
    elseif t_end < duration
        start_wf = Waveform(wf.f,duration - t_end)

        end_wf = linear_ramp(;duration=duration-t_end,
            start_value=wf(t_end),
            stop_value=end_value
        )

        return append(start_wf,end_wf)
    else
        return wf      
    end

    return new_wf

end



function main()
    wf = Waveform(t->1+t^2,1)
    max_slope = 100
    begin_value = 0.5
    end_value = 1.0
    new_wf = pin_waveform_edges(wf,max_slope,begin_value,end_value)
    println(new_wf.duration)
end


main()



