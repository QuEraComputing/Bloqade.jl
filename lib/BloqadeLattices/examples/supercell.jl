using BloqadeLattices



begin
    bounds = zeros((2,2))
    bounds[:,1] .= (1,3)
    bounds[:,2] .= (3,0)

    pp = Parallelepiped(bounds)

    x = [rand(2) for i in 1:10]

   within_cell.(pp,x)

end

