"""
obtain gradients by numerical differentiation.
"""
function numgrad(f, x::AbstractArray; δ=1e-5)
    res = similar(x)
    for i=1:length(x)
        x[i] += δ/2
        pos = f(x)
        x[i] -= δ
        neg = f(x)
        x[i] += δ/2
        res[i] = (pos-neg)/δ
    end
    res
end
