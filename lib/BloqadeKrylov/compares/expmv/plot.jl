
using Plots
using LaTeXStrings

Nt = 1000
scale = 6.0

f = open("3-by-3.$(Nt)._$(scale)_.ts","r")
ts = convert(Array{Float64,1}, collect(readeach(f,Float64)))
close(f)


f = open("3-by-3.$(Nt)._$(scale)_.exmp","r")
Ns = convert(Array{Int,1}, collect(readeach(f,Int)))
close(f)


f = open("3-by-3.$(Nt)._$(scale)_.exm_mply.ovh","r")
Nnew_ovhd = convert(Array{Int,1}, collect(readeach(f,Int)))
close(f)


f = open("3-by-3.$(Nt)._$(scale)_.exm_mply.impl","r")
Nnew_impl = convert(Array{Int,1}, collect(readeach(f,Int)))
close(f)


Plots.plot(ts, Ns,linestyle=:dash,label="expmv", ylabel = "# of mul!(::SumOfLinop,v)", xlabel=L"$t$")
Plots.plot!(ts, Nnew_ovhd,linestyle=:dash,label="expm_multiply, overhead (s,m)", ylabel = "# of mul!(::SumOfLinop,v)", xlabel=L"$t$")
Plots.plot!(ts, Nnew_impl,linestyle=:dash,label="expm_multiply, impl", ylabel = "# of mul!(::SumOfLinop,v)", xlabel=L"$t$")
Plots.plot!(ts, Nnew_impl.+Nnew_ovhd,linestyle=:dash,label="expm_multiply (total = overhead+impl)", ylabel = "# of mul!(::SumOfLinop,v)", xlabel=L"$t$")
Plots.title!("3-by-3, Nt = $Nt, scale = $scale")

savefig("3-by-3.$(Nt)._$(scale)_.png")