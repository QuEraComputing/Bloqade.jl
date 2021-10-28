using DelimitedFiles
using Sockets

"""
    cloud_f(postfunc; params_ncol, port)

* `postfunc` isa a function that take spin configurations as input, and the loss as output.
* `params_ncol` is the number of data columns.
* `port` is the socket application port.
"""
function cloud_f(postfunc; params_ncol::Int, port::Int, overwrite_file::Bool=false)
    k = 0
    function loss(x)
        c = connect(port)
        if !isopen(c)
            println("connection fail!")
        else
            println("connection success!")
        end

        # send file
        k += 1
        ifname = "Parameters_$k.dat"
        if isfile(ifname) && !overwrite_file
            error("file name conflict! got input file $ifname.")
        end
        writedlm(ifname, reshape(x,:,params_ncol))
        println(c, ifname)

        # receive file
        local ofname
        while isopen(c)
            ofname = readline(c, keep=false)
            println(stdout, ofname)
            ofname
        end
        close(c)
        println(stdout, "Get datafile: $ofname")
        if !isfile(ofname)
            error("output file not found! got output file $ofname.")
        end
        spinconfigs = readdlm(ofname, Int)
        postfunc(spinconfigs)
    end
end

function postfunc(spinconfigs)
    sum(spinconfigs)/size(spinconfigs,1)
end

loss = cloud_f(postfunc; params_ncol=3, port=2020, overwrite_file=true)
@show loss(randn(12))
@show loss(randn(12))

#using Test
#@testset "postfunc" begin
    #spins = [0 1; 1 1; 1 1]
    #@test postfunc(spins) == 5/3
#end
