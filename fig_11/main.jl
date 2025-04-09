using CSV
using Plots
using PlotlyJS
using DataFrames
using Dates
using FileIO
using ImageMagick

r = CSV.read("rateCali.csv", DataFrame) 
mj= 13

N = r.N[1] # N population size
I_scaled = N^(-1) * r.I0 #scaled infected data



plotlyjs()

p = Plots.plot(r.Fecha[(mj-10):end], r.I0[(mj-10):end],
                line=(2, :solid),
                label="Number infected per day", color="red",
                size=(675, 440), titleposition=:left,
                legend=:outerbottomleft, 
                Dates.format.(r.Fecha[mj:7:end], "yyyy-mm-dd")
                )
Plots.savefig("p.png") # fig 10
img169 = load("p.png")
save("dataset_cali.jpg",img169)