abstract type AbstractTrialState{T <: Real, C}; end

# concrete TrialStates that aren't ProductStates should have two fields
# left_flips and right_flips that are PushVector{Int} of max-length = to the number of sites
# these will cache the sites that are going to be flipped

abstract type AbstractProductState{T, C} <: AbstractTrialState{T, C}; end

function weightchange(psi::AbstractProductState{T, Bool}, initial::AbstractArray{Bool}, sites::AbstractVector{Int})::T where T
    isempty(sites) && return one(T)
    prod(s -> weightchange(psi, @inbounds initial[s]), sites)
end

function logweightchange(psi::AbstractProductState{T, Bool}, initial::AbstractArray{Bool}, sites::AbstractArray{Int})::T where T
    isempty(sites) && return zero(T)
    sum(s -> logweightchange(psi, @inbounds initial[s]), sites)
end

struct PlusState{T, C} <: AbstractProductState{T, C}; end
weightchange(psi::PlusState{T, Bool}, initial::Bool) where T = one(T)
logweightchange(psi::PlusState{T, Bool}, initial::Bool) where T = zero(T)


struct BoolDict{T} <: AbstractDict{Bool, T}
    a::Vector{T}
    BoolDict(d::AbstractDict{Bool, T}) where T = new{T}([d[false], d[true]])
end
Base.haskey(d::BoolDict, key::Bool) = true
Base.haskey(d::BoolDict, key) = false
Base.getindex(d::BoolDict{T}, key::Bool) where T = (@inbounds d.a[key + 1])
Base.get(d::BoolDict{T}, key::Bool) where T = (@inbounds d.a[key + 1])
Base.isempty(d::BoolDict) = false


struct ProductState{T, C, D <: AbstractDict{C, T}} <: AbstractProductState{T, C}
    p::D
    lp::D
    ProductState{T, C}(p::AbstractDict{C, T}) where {T, C} = new{T, C, typeof(p)}(p, Dict([k => log(v) for (k, v) in p]))
    ProductState{T, Bool}(p::AbstractDict{Bool, T}) where T = new{T, Bool, BoolDict{T}}(BoolDict(p), BoolDict(Dict([k => log(v) for (k, v) in p])))
    ProductState{T}(p::AbstractDict{C, T}) where {T, C} = ProductState{T, C}(p)
    ProductState(p::AbstractDict{C, T}) where {T, C} = ProductState{T}(p)
end
weightchange(psi::ProductState{T, C}, initial::C, final::C) where {T, C} = psi.p[final]/psi.p[initial]
weightchange(psi::ProductState{T, Bool}, initial::Bool) where T = weightchange(psi, initial, !initial)

logweightchange(psi::ProductState{T, C}, initial::C, final::C) where {T, C} = psi.lp[final] - psi.lp[initial]
logweightchange(psi::ProductState{T, Bool}, initial::Bool) where T = logweightchange(psi, initial, !initial)
