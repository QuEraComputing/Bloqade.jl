using Sockets
using DelimitedFiles

const k = Ref(1)
# a simple server (device side)
t = @async begin
    server = listen(2020)
    while true
        sock = accept(server)
        @async while isopen(sock)
            s = readline(sock, keep=false)
            if isempty(s)
                continue
            end
            println(stdout, "Get file $s...")
            if !isfile(s)
                println(stdout, "$s is not a file!")
                write(sock, "not a file!")
                continue
            end
            params = readdlm(s)
            println(stdout, "Get params $params, running the experiment...")
            # 20 spins, 200 measurements
            result = rand([-1,1], 200, 20)
            ofname = "Measurementdata_$(k[]).dat"
            writedlm(ofname, result)
            write(sock, ofname)
            println(stdout, "Data sent!")
            close(sock)
            k[] += 1
        end
    end
end

wait(t)
