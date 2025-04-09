using CSV
using Plots
using PlotlyJS
using DataFrames
using Dates
using FileIO
using ImageMagick

include("./iran_case.jl")
include("./iran_case2.jl")



############################################################################
############################################################################

################## Rolling historical data to one forecast 

#size_forecast = 7 # testing
size_day = 2 # testing
start_forecast = 11 # testing
size_rolling_historical_window = 10
gamma = 1 / 14
#################### time data subdivitions to use Euler Maruyama
replications = 1000 #testing
subdivitions = 1000 # 50 # number of subdivitions into step values
mean_noise = 0.0 # testing
alp = 0.05 

################# Preparing data
#r = CSV.read("rateIran_paper.csv", DataFrame) #           
r = CSV.read("rateCali.csv", DataFrame) # forecast paper until 15 sept - change Fecha
 
N = r.N[1] # N population size
I_scaled = N^(-1) * r.I0 #scaled infected data

################## 
extrem_forecast = zeros((length(I_scaled)-size_rolling_historical_window - 1), 2)
mean_forecast = zeros((length(I_scaled)-size_rolling_historical_window - 1), 1)
rmse_forecast = zeros((length(I_scaled)-size_rolling_historical_window - 1), 1)
end_forecast = (length(I_scaled)-size_rolling_historical_window)
for j in 1:(end_forecast -1) 
        model = Iran_get_forecast(start_forecast + (j - 1) - size_rolling_historical_window ,start_forecast + (j - 1),2,size_day,replications,subdivitions,alp,mean_noise,gamma)
        Trajectories = model[1]
        mean_trajectories = model[2]
        rmse_forecast[j] = model[3]
 
            dat = Trajectories[:, 2]
            mean_forecast[j, 1] = mean(dat) #quantile(dat, 0.5)#0.95
            extrem_forecast[j, 1] = minimum(dat)
            extrem_forecast[j, 2] = maximum(dat)
        
end

r = CSV.read("rateCali.csv", DataFrame) # forecast paper until 15 sept - change Fecha
 
#################### eliminate null Infected values
Float_I = Float64.(r.I0) 
####################
# scaling data

plotlyjs()
# xticks=(r.Fecha, Dates.format.(r.Fecha, "yyyy-mm-dd"))
g = Plots.plot(r.Fecha, Float_I,
    line=(2, :solid),
    label="Daily infected", color="red",
    size=(1495, 640), titleposition=:left,
    legend=:outerbottomleft,
    xticks=(r.Fecha[1:15:end], 
            Dates.format.(r.Fecha[1:15:end], "yyyy-mm-dd")),
    xrotation = 90,
    xtickfont=(8),
    widen=false,
    tickdirection=:out)



g = Plots.plot!(g,r.Fecha[(start_forecast + 1):end],mean_forecast*N,
    line=(2, :solid),
    label="Mean forecast Itô-SIS ", color="blue",
    size=(1495, 640), titleposition=:left)
    #p = Plots.plot!(p,r.Fecha[j:(j+size_forecast-1)],mean_trajectories[:, 1]*N,
    #line=(2, :solid),
    #label="Mean forecast Itô-SIS", color="orange",
    #size=(675, 200), titleposition=:left,
    #legend=:bottomright)
g = Plots.plot(g,r.Fecha[(start_forecast+1):end], extrem_forecast[:, 1]*N,
    fillrange=extrem_forecast[:, 2]*N, fillalpha=0.15, c=1,
    label="Forecast Itô-SIS band  ")
g  
Plots.savefig("c_rolling.png") # fig 6
img26 = load("c_rolling.png")
save("c_rolling.jpg",img26)











