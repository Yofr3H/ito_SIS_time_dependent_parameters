# Gonorrhea proofs
include("./gn_case.jl")
using CSV
using Plots
using PlotlyJS
using DataFrames
using Dates
using FileIO
using ImageMagick


r = CSV.read("dat_gn.csv", DataFrame) # forecast paper until 15 sept - change Fecha

    
    plotlyjs()
    # xticks=(r.Fecha, Dates.format.(r.Fecha, "yyyy-mm-dd"))
    p = Plots.plot(r.Year, r.I,
        line=(2, :solid),
        label="Number infected per day", color="red",
        size=(675, 200), titleposition=:left,
        legend=:topright)
    # scaling data

#I_scaled = r.I



size_forecast = 10 # testing
size_day = 2 # testing
mj= 50 # testing
j = 60 # testing
gamma = 1/55
#################### time data subdivitions to use Euler Maruyama
replications = 1000 #testing
subdivitions = 1000 # 50 # number of subdivitions into step values
mean_noise = 0.0 # testing
alp = 0.05 # percentil to cut outiliers


#M = SDE_SISg(subdivitions, size_forecast, size_day, mj, j, mean_noise, gamma, alp, I_scaled)
        

model = model = gn_get_forecast(mj,j,size_forecast,size_day,replications,subdivitions,alp,mean_noise,gamma) # quantile procedure

Trajectories = model[1]
mean_trajectories = model[2]
RMSE_mean_estimated = model[3]


################## plot forecast and comparision MLP method j-10, j= 178


Time = range(j, j + size_forecast - 1, size_forecast)
extrem_trajectories = zeros(length(Time), 2)
median_trajectories = zeros(length(Time))
for i = 1:(length(Time))
    dat = Trajectories[:, i]
    median_trajectories[i, 1] = quantile(dat, 0.5)#0.95
    extrem_trajectories[i, 1] = minimum(dat)
    extrem_trajectories[i, 2] = maximum(dat)
end




plotlyjs()

g1 = Plots.plot(r.Year[(mj-10):end], r.I_1[(mj-10):end],
    line=(2, :solid),
    label="Number infected per day", color="red",
    size=(675, 440), titleposition=:left,
    legend=:outerbottomleft,
    xrotation = 90,
    xtickfont=(8),
    widen=false,
    tickdirection=:out)
     #xticks=(range(r.Fecha[mj],r.Fecha[end],8), "yy-mm-dd"))
N = mean(r.Nt)

g1 = Plots.plot!(g1,r.Year[j:(j+size_forecast-1)],median_trajectories[:, 1]*N,
    line=(2, :solid),
    label="Median forecast Itô-SIS", color="blue",
    size=(675, 440), titleposition=:left)
    #p = Plots.plot!(p,r.Fecha[j:(j+size_forecast-1)],mean_trajectories[:, 1]*N,
    #line=(2, :solid),
    #label="Mean forecast Itô-SIS", color="orange",
    #size=(675, 200), titleposition=:left,
    #legend=:bottomright)
g1 = Plots.plot(g1,r.Year[j:(j+size_forecast-1)], extrem_trajectories[:, 1]*N,
    fillrange=extrem_trajectories[:, 2]*N, fillalpha=0.15, c=1,
    label="forecast Itô-SIS band")
g1
Plots.savefig("g1.png")# pronostico con datos entre 1990 y 2000 para obtener 10 años de pronostico
img213 = load("g1.png")
save("g1.jpg",img213)

#testing beta behaviour

scaled_beta = beta.^(-1) * 365.25^(-1)
g2 = Plots.plot(r.Year[2:80],scaled_beta[2:80]) 
