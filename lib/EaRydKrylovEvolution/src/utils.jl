function unsafe_log2i(x)
    error("unsupported type $(typeof(x))")
end

for N in [8, 16, 32, 64, 128]
    T = Symbol(:Int, N)
    @eval begin
        unsafe_log2i(x::$T) = $(N - 1) - leading_zeros(x)
    end
end
