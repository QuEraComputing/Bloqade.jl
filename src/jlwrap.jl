using PythonCall: pysysmodule, pyosmodule

const pybloqade_module = PythonCall.pynew()

function init_bloqade()
    jl = pybloqade_module
    sys = pysysmodule
    os = pyosmodule
    if C.CTX.is_embedded
        # in this case, Julia is being embedded into Python by juliacall, which already exists
        pycopy!(jl, sys.modules["bloqade"])
    end
end
