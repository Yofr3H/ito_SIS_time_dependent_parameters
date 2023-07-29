using CSV
using Plots
using Plots
using PlotlyJS
using DataFrames

include("./iran_forward_case.jl")
include("./iran_backward_case.jl")

################## Generating forward values Table 1

size_historical_data = 9  # number of historical data to use: 8,9,10,29,30,31,60
I_k = 178 # (178 - date: Agust 15, 2020)
size_forecast = 30 # size_forecast
size_rolling_window = 0 

forward_case = Iran_forward(size_historical_data,I_k,size_forecast,size_rolling_window)
Plot_original_data = forward_case[3]
RMSE_mean_estimated = forward_case[1] # approximation: 240.7

#Table 1
#size_historical_data   RMSE_mean_estimated 
#8                      255.7446
#9                      240.6787
#10                     306.2132
#29                     428.4117
#30                     407.1304
#31                     432.6823
#60                     1528.5964

################## plot forecast and comparision MLP method

Trajectories = forward_case[2]
Time = range(j, j + size_forecast - 1, size_forecast)
extrem_trajectories = zeros(length(Time), 2)
median_trajectories = zeros(length(Time))
for i = 1:(length(Time))
    dat = Trajectories[:, i]
    median_trajectories[i, 1] = quantile(dat, 0.5)#0.95
    extrem_trajectories[i, 1] = minimum(dat)
    extrem_trajectories[i, 2] = maximum(dat)
end

d = CSV.read("dat_paper.csv", DataFrame)
aux2 = RSME(d.fore,I_scaled[(length(I_scaled) - size_forecast ):(end)]*N) # 407.60

r = CSV.read("rateIran_paper.csv", DataFrame) # forecast paper until 15 sept - change Fecha
plotlyjs()
p = Plots.plot(r.Fecha, r.I0,
    line=(2, :solid),
    label="Number infected per day", color="red",
    size=(675, 200), titleposition=:left,
    legend=:topright)

N = r.N[1] # N population size
I_scaled = N^(-1) * r.I0 #scaled infected data
p = Plots.plot!(p,r.Fecha[178:208],d.fore,
    line=(2, :solid),
    label="MLP forecast", color="green",
    size=(675, 200), titleposition=:left,
    legend=:bottomright) #case iran only

p = Plots.plot!(p,r.Fecha[j:(j+size_forecast-1)],median_trajectories[:, 1]*N,
    line=(2, :solid),
    label="Median forecast Itô-SIS", color="blue",
    size=(675, 200), titleposition=:left,
    legend=:bottomright)
p = Plots.plot(p,r.Fecha[j:(j+size_forecast-1)], extrem_trajectories[:, 1]*N,
    fillrange=extrem_trajectories[:, 2]*N, fillalpha=0.15, c=1,
    label="forecast Itô-SIS band", legend=:topleft)
p

################ Boxplot of RMSE median, using the rolling window

size_rolling_window = 30
data_boxplot_forward_case = Iran_forward(mj,j,size_forecast,size_rolling_window)

data_boxplot_forward = filter(!iszero, data_boxplot_forward_case[4]*N)
trace1 = box(
    x=data_boxplot_forward,
    boxpoints="all",
    name="RMSE"
    )
PlotlyJS.plot([trace1]) # conicides with boxplot Figure 8


################## Generating backward values Table 1

size_historical_data = 8  # number of historical data to use: 8,14,15,16,30,60
I_k = 178 # (178 - Agust 15, 2020)
size_forecast = 30 # size_forecast
size_rolling_window = 0 
backward_case = Iran_backward(size_historical_data,I_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = backward_case[1] # approximation: 414.1 

#Table 1 (continuation)
#size_historical_data   RMSE_mean_estimated 
#8                      414.1980
#14                     419.7058
#15                     374.3958
#16                     370.6485
#30                     396.5671
#60                     1525.0232
