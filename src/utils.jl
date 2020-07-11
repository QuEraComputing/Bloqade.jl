function unsafe_log2i end

for N in [8, 16, 32, 64, 128]
    T = Symbol(:Int, N)
    UT = Symbol(:UInt, N)
    @eval begin
        unsafe_log2i(x::$T) = $(N - 1) - leading_zeros(x)
        unsafe_log2i(x::$UT) = $(N - 1) - leading_zeros(x)
    end
end
