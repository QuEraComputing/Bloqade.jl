"""
Aquila noise model parameters in JSON format
"""
AQUILA = """{
    "coherent":{
        "local":{
            "site x":0.05,"site y":0.05,"detuning":0.5,"Rabi frequency":0.018
        },
        "global":{
            "detuning":0.18,"Rabi frequency":0.008
        }
    },
    "incoherent":{
        "dephasing rate":0.02,"relaxation rate":0.01
    },
    "readout":{
        "single atom":{"p01":0.01,"p10":0.08}
    }
}"""

"""
    function load_error_model

Load a generic Rydberg atom error model from a JSON string

# Arguments
- `config`: JSON string with error model parameters (example in src/noise_models.jl)
"""
function load_error_model(props)
    relax_op = (X+im*Y)/2 #|0⟩⟨1|

    p01 = props["readout"]["single atom"]["p01"]
    p10 = props["readout"]["single atom"]["p10"]
    relaxation_rate = props["incoherent"]["relaxation rate"]
    dephasing_rate = props["incoherent"]["dephasing rate"]
    δΩrel = props["coherent"]["global"]["Rabi frequency"]
    δΔ = props["coherent"]["global"]["detuning"]
    δx = props["coherent"]["local"]["site x"]
    δy = props["coherent"]["local"]["site y"]
    δΔ_inhom = props["coherent"]["local"]["detuning"]
    δΩ_inhom = props["coherent"]["local"]["Rabi frequency"] #changed from .02 (number density is sensitive to this parameter)

    function coherent_noisy(h)
        (atoms,ϕ,Ω,Δ) = get_rydberg_params(h)
        if Δ === nothing
            throw("Δ needs to be specified!")
        end
        if eltype(atoms) == Tuple{Float64}
            atoms = [(first(a), 0.0) for a in atoms]
        end

        function sample_noisy_ham() #return a function that modifies h
            randomize((x,y)) = (x+δx * randn(), y + δy*randn())
            atoms_noisy = randomize.(atoms) #randomize atom positions
            #add coherent drift in Ω and inhomogeneity
            Ω_noisy = Ω .* (1+δΩrel*randn() .+ δΩ_inhom * randn(length(atoms)))
            #add coherent drift in Δ and inhomogeneity
            Δ_noisy = δΔ*randn()+Δ .+ δΔ_inhom * randn(length(atoms))
            return rydberg_h(
                atoms_noisy;
                Ω = Ω_noisy,
                Δ = Δ_noisy,
                ϕ = ϕ
            )
        end
    end

    function collapse_operators(nqubits)
        [[
            sqrt(dephasing_rate) * SparseMatrixCSC(mat(put(nqubits, q => Z))) #single-qubit relaxation
            for q in 1:nqubits
        ];
        [
            sqrt(relaxation_rate) * mat(put(nqubits, q => relax_op)) #single-qubit relaxation
            for q in 1:nqubits
        ]]
    end

    function confusion_mat(N)
        M = [[1-p01 p10]; [p01 1-p10]]
        return kronecker([M for i in 1:N]...) #readout
    end

    ErrorModel(
        confusion_mat,
        collapse_operators,
        coherent_noisy 
    )
end

#add waveforms and numbers (required for Hamiltonian sampling)
function Base.:+(a::Float64, b::Waveform)
    Waveform(b.duration) do t
        b(t)+a
    end
end

function Base.:+(a::Waveform, b::Float64)
    Waveform(a.duration) do t
        a(t)+b
    end
end


"""
function Aquila

    Create an ErrorModel representing the noise model of Aquila.
"""
Aquila() = load_error_model(JSON.parse(AQUILA))