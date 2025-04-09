using CSV
using Plots
using PlotlyJS
using DataFrames
using Dates
using FileIO
using ImageMagick


r = CSV.read("rateIran_paper.csv", DataFrame)
mj= 78

N = r.N[1] # N population size
I_scaled = N^(-1) * r.I0 #scaled infected data

p = Plots.plot(r.Fecha[(mj-10):end], r.I0[(mj-10):end],
    line=(2, :solid),
    label="Number infected per day", color="red",
    size=(675, 440), titleposition=:left,
    legend=:outerbottomleft,
    xticks=(r.Fecha[mj:7:end], 
            Dates.format.(r.Fecha[mj:7:end], "yyyy-mm-dd")),
    xrotation = 90,
    xtickfont=(8),
    widen=false,
    tickdirection=:out)

d = CSV.read("dat_paper.csv", DataFrame)

p = Plots.plot!(p,r.Fecha[178:208],d.fore,
    line=(2, :solid),
    label="MLP forecast", color="green",
    size=(675, 440), titleposition=:left) #case iran only

Plots.savefig("iran.png")
img26 = load("iran.png")
save("iran.jpg",img26)

