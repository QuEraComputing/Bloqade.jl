function emulate_routine!(r::RydbergReg{ComplexLayout}, duration::Number, h::AbstractTerm, cache::AbstractMatrix)
    st = vec(r.state)
    update_term!(cache, h, r.subspace)
    expmv!(-im * duration, cache, st)
    return r
end

function emulate_routine!(r::RydbergReg{RealLayout}, duration::Number, h::AbstractTerm, cache::AbstractMatrix)
    error("discrete emulator does not support RealLayout yet")
end

function emulate_routine!(r::Yao.ArrayReg, duration::Number, h::AbstractTerm, cache::AbstractMatrix)
    st = vec(r.state)
    update_term!(cache, h)
    expmv!(-im * duration, cache, st)
    return r
end
