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
struct CFET2_1 <: CFETTables
    xs::Vector{Float64}
    Gs::Vector{Vector{Float64}}
    function CFET2_1()
        xs =  [0.5]
        Gs = [[1.0]]
        return new(xs,Gs)
    end
end

# CFET with 4th order with 2 exponentials:
struct CFET4_2 <: CFETTables
    xs::Vector{Float64}
    Gs::Vector{Vector{Float64}}
    function CFET4_2()
        xs =  [1/2 - √3/6, 1/2 + √3/6]
        Gs = [[1/4 + √3/6, 1/4 - √3/6], 
              [1/4 - √3/6, 1/4 + √3/6]]
        return new(xs,Gs)
    end
end

# CFET with 4th order with 2 exponentials:
#= tested! auto matic generate way.
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
struct CFET6_5  <: CFETTables
    xs::Vector{Float64}
    Gs::Vector{Vector{Float64}}
    function CFET6_5()
        
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


# CFET with 8th order with 11 exponentials:
# (unoptimized coefficients)
struct CFET8_11  <: CFETTables
    xs::Vector{Float64}
    Gs::Vector{Vector{Float64}}
    function CFET8_11()
        
        xs, ws = GaussQuadrature.legendre(4) # this can change, but generally choose s = m is ok
        xs .+= 1.0
        xs ./= 2.0
        ws ./= 2.0 

        # Fs is the coefficients matching convention from Ref[1]
        # the second index runs from 1 ~ N/2 
        f11 = 0.169715531043933180094151;  f12 = 0.152866146944615909929839
        f13 = 0.119167378745981369601216;  f14 = 0.068619226448029559107538
        f21 = 0.379420807516005431504230;  f22 = 0.148839980923180990943008
        f23 = −0.115880829186628075021088; f24 = −0.188555246668412628269760
        f31 = 0.469459306644050573017994;  f32 = −0.379844237839363505173921
        f33 = 0.022898814729462898505141;  f34 = 0.571855043580130805495594
        f41 = −0.448225927391070886302766; f42 = 0.362889857410989942809900
        f43 = −0.022565582830528472333301; f44 = −0.544507517141613383517695
        f51 = −0.293924473106317605373923; f52 = −0.026255628265819381983204
        f53 = 0.096761509131620390100068;  f54 = 0.000018330145571671744069
        f61 = 0.447109510586798614120629;  f63 = −0.200762581179816221704073

        Fs = [[f11, f12, f13, f14], 
              [f21, f22, f23, f24],
              [f31, f32, f33, f34],
              [f41, f42, f43, f44],
              [f51, f52, f53, f54],
              [f61,   0, f63,   0]] # the zeros here are requred by reveral of T-symm

        
        # make the rest base on T-symm condition (R1)
        for p in 5:-1:1
            push!(Fs,[ (n+1)%2==0 ? Fs[p][n] : -Fs[p][n] for n in 1:length(Fs[p])])
        end 


        return new(xs,__get_Gs(Fs, xs, ws))
    end
end

