
using Bloqade
using BloqadeKrylov
using Plots
using LaTeXStrings

f = open("Ns","r")
Ns = convert(Array{Float64,1}, collect(readeach(f,Int64)))
close(f)

f = open("es_kry","r")
reskry = convert(Array{Float64,1}, collect(readeach(f,Float64)))
close(f)

f = open("es_m42","r")
res42 = convert(Array{Float64,1}, collect(readeach(f,Float64)))
close(f)

f = open("es_m65","r")
res65 = convert(Array{Float64,1}, collect(readeach(f,Float64)))
close(f)

f = open("es_m811","r")
res811 = convert(Array{Float64,1}, collect(readeach(f,Float64)))
close(f)

Plots.plot(Ns, reskry, marker=".",label="krylov", xaxis=:log, yaxis=:log, xlabel = L"$N_T$", ylabel=L"$ϵ$")
Plots.plot!(Ns, res42, marker=".",label="magnus42", xaxis=:log, yaxis=:log, xlabel = L"$N_T$", ylabel=L"$ϵ$")
Plots.plot!(Ns, res65, marker=".",label="magnus65", xaxis=:log, yaxis=:log, xlabel = L"$N_T$", ylabel=L"$ϵ$")
Plots.plot!(Ns, res811, marker=".",label="magnus811", xaxis=:log, yaxis=:log, xlabel = L"$N_T$", ylabel=L"$ϵ$")
Plots.title!("3x3 Square")
savefig("adb.png")

Plots.plot(Ns, reskry, marker=".",label="krylov", xaxis=:log, yaxis=:log, xlabel = L"$N_T \times s (\# of expmv)$", ylabel=L"$ϵ$")
Plots.plot!(Ns.*2, res42, marker=".",label="magnus42", xaxis=:log, yaxis=:log, xlabel = L"$N_T \times s (\# of expmv)$", ylabel=L"$ϵ$")
Plots.plot!(Ns.*5, res65, marker=".",label="magnus65", xaxis=:log, yaxis=:log, xlabel = L"$N_T \times s (\# of expmv) $", ylabel=L"$ϵ$")
Plots.plot!(Ns.*11, res811, marker=".",label="magnus811", xaxis=:log, yaxis=:log, xlabel = L"$N_T \times s (\# of expmv)$", ylabel=L"$ϵ$")
Plots.title!("3x3 Square")
savefig("adb_scale.png")