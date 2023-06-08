# arxiv:1102.5071 eq.(59)
# notation: 
#   Gs: number of exponential time propogator (stored in the order from right to left)
#        Gs has size (# of exponential time propogator,length(xs))
#   xs: location of the gaussian lagendre quadrature
struct CFET42 <: CFETTables
    xs::Vector{Float64}
    Gs::Vector{Vector{Float64}}
    function CFET42()
        xs =  [1/2 - √3/6, 1/2 + √3/6]
        Gs = [[1/4 + √3/6, 1/4 - √3/6], 
              [1/4 - √3/6, 1/4 + √3/6]]
        return new(xs,Gs)
    end
end