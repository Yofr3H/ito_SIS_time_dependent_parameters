using CSV
using Plots
using PlotlyJS
using DataFrames
using Dates
using FileIO
using ImageMagick

include("./iran_case.jl")
#include("./iran_case2.jl")



############################################################################
############################################################################

################## Rolling historical data to one forecast 


size_day = 5 # testing
start_forecast = 160 # testing 31: 2020-04-01 cali
size_rolling_historical_window = 30
gamma = 1 / 14
#################### time data subdivitions to use Euler Maruyama
replications = 1000 #testing
subdivitions = 1000 # 50 # number of subdivitions into step values
mean_noise = 0.0 # testing
alp = 0.05 

################# Preparing data
#r = CSV.read("rateIran_paper.csv", DataFrame) #           
#r = CSV.read("rateCali.csv", DataFrame) # forecast paper until 15 sept - change Fecha
 
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
# forecast paper until 15 sept - change Fecha
 
#################### eliminate null Infected values
Float_I = Float64.(r.I0) 
####################
# scaling data

plotlyjs()


#plot rolling between specific dates
date_first = Date("2020-08-05")  # example cali: 2020-05-31 
date_end = Date("2020-09-05") # example cali: 2020-07-15

idx_first = findfirst(==(date_first), r.Fecha)
idx_end = findfirst(==(date_end), r.Fecha)

g = Plots.plot(r.Fecha[idx_first:idx_end], Float_I[idx_first:idx_end],
    line=(2, :solid),
    label="Daily infected", color="red",
    size=(1495, 640), titleposition=:left,
    legend=:outerbottomleft,
    xticks=(r.Fecha[idx_first:5:idx_end], 
            Dates.format.(r.Fecha[idx_first:5:idx_end], "yyyy-mm-dd")),
    xrotation=90,
    xtickfont=(8),
    widen=false,
    tickdirection=:out)

range_forecast = start_forecast:(start_forecast + size(mean_forecast, 1) - 1)
forecast_dates = r.Fecha[range_forecast]
forecast_mask = (forecast_dates .>= date_first) .& (forecast_dates .<= date_end)
valid_indices = findall(forecast_mask)

g = Plots.plot!(g, forecast_dates[valid_indices], mean_forecast[valid_indices] * N,
    line=(2, :solid),
    label="Mean forecast Itô-SIS", color="blue",
    size=(1495, 640), titleposition=:left)

g = Plots.plot!(g, forecast_dates[valid_indices],
    extrem_forecast[valid_indices, 1] * N,
    fillrange=extrem_forecast[valid_indices, 2] * N,
    fillalpha=0.15, c=1,
    label="Forecast Itô-SIS band")
g  
Plots.savefig("c_rolling.png") # fig 6
img26 = load("c_rolling.png")
save("c_rolling_Milstein_1.jpg",img26)

