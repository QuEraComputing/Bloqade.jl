# Ref[1] arxiv:1102.5071 eq.(59)
# notation: 
#   
#   Ωi = ∑_(m=1) Gs[i,m] Ham(t + xs[m]δt)   
#
#   Gs: number of exponential time propogator (stored in the order from right to left)
#        Gs has size (# of exponential time propogator,length(xs))
#   xs: location of the gaussian lagendre quadrature (which we choosen to be N/2, but generally can be larger)
#       

# CFET with 2nd order with 1 exponentials:
# (2nd order midpoint rule)
struct CFET21 <: CFETTables
    xs::Vector{Float64}
    Gs::Vector{Vector{Float64}}
    function CFET42()
        xs =  [0.5]
        Gs = [[1.0]]
        return new(xs,Gs)
    end
end

# CFET with 4th order with 2 exponentials:
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

# CFET with 4th order with 2 exponentials:
#= tested!
struct CFET42Gen <: CFETTables
    xs::Vector{Float64}
    Gs::Vector{Vector{Float64}}
    function CFET42Gen()
        xs, ws = GaussQuadrature.legendre(2) 
        xs .+= 1.0
        xs ./= 2.0
        ws ./= 2.0

        Fs = [[1/2, 1/3 ],
              [1/2, -1/3]]
        
        Gs = __get_Gs(Fs,xs,ws) 
        return new(xs,Gs)
    end
end
=#


# CFET with 6th order with 5 exponentials:
# (unoptimized coefficients)
struct CFET65  <: CFETTables
    xs::Vector{Float64}
    Gs::Vector{Vector{Float64}}
    function CFET65()
        
        xs, ws = GaussQuadrature.legendre(5) # this can change, but generally choose s = m is ok
        xs .+= 1.0
        xs ./= 2.0
        ws ./= 2.0 

        # Fs is the coefficients matching convention from Ref[1]
        # the second index runs from 1 ~ N/2 
        Fs = [[0.16                  , 0.14587456942714338561, 0.11762370828143015682], 
              [0.38752405202531186588, 0.15089113704380764664, −0.12805075909013044594]]
        
        
        push!(Fs,[1-2*Fs[2][1]-2*Fs[1][1], 0 , -2*Fs[2][3]-2*Fs[1][3]])
        
        # make the rest base on T-symm condition (R1)
        push!(Fs,[ (n+1)%2==0 ? Fs[2][n] : -Fs[2][n] for n in 1:length(Fs[2])])
        push!(Fs,[ (n+1)%2==0 ? Fs[1][n] : -Fs[1][n] for n in 1:length(Fs[1])])

        return new(xs,__get_Gs(Fs, xs, ws))
    end
end


