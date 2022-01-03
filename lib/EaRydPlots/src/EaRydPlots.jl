module EaRydPlots

using Makie
using EaRydCore

export bitstring_histgram, bitstring_histgram!

bitstring_histgram() = Makie.not_implemented_for(bitstring_histgram)
const BitStringHistgram{ArgType} = Makie.Combined{bitstring_histgram, ArgType}

Makie.plotsym(::Type{<:BitStringHistgram}) = begin
    :BitStringHistgram
end

function bitstring_histgram(args...; attributes...)
    Makie.plot(BitStringHistgram, args...; attributes...)
end

function bitstring_histgram!(args...; attributes...)
    Makie.plot!(BitStringHistgram, args...; attributes...)
end

function Makie.default_theme(scene, ::Makie.Type{<:BitStringHistgram})
    Attributes(nlargest_possible=10)
end

function find_largest(probs, num)
    if length(probs) ≤ num
        return 1:length(probs)
    end

    perm = sortperm(probs, rev=true)
    indices = perm[1:num]
    return indices
end

function Makie.plot!(bh::BitStringHistgram{<:Tuple{<:ArrayReg}})
    reg = to_value(bh[1])
    nlargest_possible = to_value(bh.nlargest_possible)
    state = statevec(reg)
    probs = abs2.(state)
    indices = find_largest(probs, nlargest_possible)
    barplot!(bh, 1:length(indices), probs[indices], bar_labels=indices .- 1)
    return bh
end

function Makie.plot!(bh::BitStringHistgram{<:Tuple{<:RydbergReg}})
    reg = to_value(bh[1])
    nlargest_possible = to_value(bh.nlargest_possible)
    state = statevec(reg)
    probs = abs2.(state)
    indices = find_largest(probs, nlargest_possible)

    bitstrings = reg.subspace.subspace_v[indices]
    barplot!(bh, 1:length(bitstring), probs[indices]; bar_labels=bitstrings)
    return bh
end

end # EaRydPlots
