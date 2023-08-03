using CSV
using Plots
using Plots
using PlotlyJS
using DataFrames

include("./delhi_adj_forward_case.jl")
include("./delhi_adj_backward_case.jl")

# ... PLEASE ... path must be modified

################## Generating forward values Table 1

size_historical_data = 6   # number of historical data to use: 8,9,10,29,30,31,60
t_k = 44 ## (t_k > size_historical_data + 1) (44 - date: April 26, 2020)
size_forecast = 20 #20 #30 # size_forecast
size_rolling_window = 0 

forward_case = Delhi_adj_forward(size_historical_data,t_k,size_forecast,size_rolling_window)
Plot_original_data = forward_case[3]
RMSE_mean_estimated = forward_case[1] # approximation: 1182.1214

#Delhi case - adjust data 18 and 19 April 2020 (size_forecast = 20)
#size_historical_data   RMSE_mean forward
#6                      1182.1214
#7                      1229.6484





################## plot forecast forward and comparision MLP method

Trajectories = forward_case[2]
j = t_k
Time = range(j, j + size_forecast - 1, size_forecast)
extrem_trajectories = zeros(length(Time), 2)
median_trajectories = zeros(length(Time))
for i = 1:(length(Time))
    dat = Trajectories[:, i]
    median_trajectories[i, 1] = quantile(dat, 0.5)#0.95
    extrem_trajectories[i, 1] = minimum(dat)
    extrem_trajectories[i, 2] = maximum(dat)
end

r= CSV.read("D:\\2021\\W_flor\\Euler_Maruyama\\SDE_SIS_13ene\\Iran_case\\rate_delhi_adj.csv", DataFrame)
    
N = r.N[1] # N population size
I_scaled = N^(-1) * r.I0 #scaled infected data

plotlyjs()
p = Plots.plot(r.Fecha[150:end], r.I0[150:end],
    line=(2, :solid),
    label="Number infected per day", color="red",
    size=(335, 200), titleposition=:left,
    legend=:topright)



p = Plots.plot!(p,r.Fecha[j:(j+size_forecast-1)],median_trajectories[:, 1]*N,
    line=(2, :solid),
    label="Median forecast Itô-SIS", color="blue",
    size=(675, 200), titleposition=:left,
    legend=:bottomright)
p = Plots.plot(p,r.Fecha[j:(j+size_forecast-1)], extrem_trajectories[:, 1]*N,
    fillrange=extrem_trajectories[:, 2]*N, fillalpha=0.15, c=1,
    label="forecast Itô-SIS band", legend=:topleft)
p


################## Generating backward values Table 1

size_historical_data = 17 # number of historical data to use: 8,14,15,16,30,60
t_k = 44 # (44 - April 26, 2020)
size_forecast = 20 # size_forecast
size_rolling_window = 0 
backward_case = Delhi_adj_backward(size_historical_data,t_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = backward_case[1] # approximation: 244.5951

#Table 1 (continuation)
#size_historical_data   RMSE_mean backward 
#8                      5095.2528
#17                     19704.1108


################## plot forecast forward

Trajectories = backward_case[2]
j = t_k
Time = range(j, j + size_forecast - 1, size_forecast)
extrem_trajectories = zeros(length(Time), 2)
median_trajectories = zeros(length(Time))
for i = 1:(length(Time))
    dat = Trajectories[:, i]
    median_trajectories[i, 1] = quantile(dat, 0.5)#0.95
    extrem_trajectories[i, 1] = minimum(dat)
    extrem_trajectories[i, 2] = maximum(dat)
end


N = r.N[1] # N population size
I_scaled = N^(-1) * r.I0 #scaled infected data

p = Plots.plot(r.Fecha[170:end], r.I0[170:end],
    line=(2, :solid),
    label="Number infected per day", color="red",
    size=(400, 190), titleposition=:right,
    legend=:outertopright)


p = Plots.plot!(p,r.Fecha[j:(j+size_forecast-1)],median_trajectories[:, 1]*N,
    line=(2, :solid),
    label="Median forecast Itô-SIS", color="blue",
    size=(400, 190),  titleposition=:right,
    legend=:outertopright)
p = Plots.plot(p,r.Fecha[j:(j+size_forecast-1)], extrem_trajectories[:, 1]*N,
    fillrange=extrem_trajectories[:, 2]*N, fillalpha=0.15, c=1,size=(400, 190),
    label="forecast Itô-SIS band",titleposition=:right,legendfontsize=7,
    legend=:none)
p

################ Boxplot of RMSE median, using the rolling window - forward case
size_historical_data = 7 # number of historical data to use: 8,14,15,16,30,60
t_k = 44 # (44 - April 26, 2020)
size_forecast = 20 # size_forecast
size_rolling_window = 7 #
data_boxplot_forward_case = Delhi_adj_forward(size_historical_data,t_k,size_forecast,size_rolling_window)

data_boxplot_forward = filter(!iszero, data_boxplot_forward_case[4]*N)
trace1 = box(
    x=data_boxplot_forward,
    boxpoints="all",
    name="RMSE"
    )
PlotlyJS.plot([trace1])  

################ Boxplot of RMSE median, using the rolling window - backward case

size_historical_data = 7 # number of historical data to use: 8,14,15,16,30,60
t_k = 44 # (44 - April 26, 2020)
size_forecast = 20 # size_forecast
size_rolling_window = 7 #
data_boxplot_backward_case = Delhi_adj_backward(size_historical_data,t_k,size_forecast,size_rolling_window)

data_boxplot_backward = filter(!iszero, data_boxplot_backward_case[4]*N)
trace1 = box(
    x=data_boxplot_backward,
    boxpoints="all",
    name="RMSE"
    )
PlotlyJS.plot([trace1]) 