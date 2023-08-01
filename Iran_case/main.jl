using CSV
using Plots
using Plots
using PlotlyJS
using DataFrames

include("./iran_forward_case.jl")
include("./iran_backward_case.jl")

################## Generating forward values Table 1

size_historical_data = 12 #7  # number of historical data to use: 8,9,10,29,30,31,60
t_k = 178 # (t_k > size_historical_data + 1) 178 (178 - date: Agust 15, 2020)
size_forecast = 30 #30 # size_forecast
size_rolling_window = 0 

forward_case = Iran_forward(size_historical_data,t_k,size_forecast,size_rolling_window)
Plot_original_data = forward_case[3]
RMSE_mean_estimated = forward_case[1] # approximation: 386.13356

#Table 1
#size_historical_data   RMSE_mean_estimated 
#6                                (1041.7964)
#7                                (1158.5953)
#8                      255.7446  (1231.74026)
#9                      240.6787  (1012.0450)
#10                     306.2132   (860.65037)
#11                                (595.1842)
#12                                (386.13356)
#13                                (1152.7663)
#29                     428.4117    (570.2873)
#30                     407.1304    (713.4612)
#31                     432.6823    (454.8572)
#60                     1528.5964   (1470.0902)

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

r = CSV.read("rateIran_paper.csv", DataFrame) # forecast paper until 15 sept - change Fecha
N = r.N[1] # N population size
I_scaled = N^(-1) * r.I0 #scaled infected data
d = CSV.read("dat_paper.csv", DataFrame)
aux2 = RSME(d.fore,I_scaled[(length(I_scaled) - size_forecast ):(end)]*N) # 407.60


plotlyjs()
p = Plots.plot(r.Fecha, r.I0,
    line=(2, :solid),
    label="Number infected per day", color="red",
    size=(675, 200), titleposition=:left,
    legend=:topright)


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

size_rolling_window = 14 #30
data_boxplot_forward_case = Iran_forward(size_historical_data,I_k,size_forecast,size_rolling_window)

data_boxplot_forward = filter(!iszero, data_boxplot_forward_case[4]*N)
trace1 = box(
    x=data_boxplot_forward,
    boxpoints="suspectedoutliers",
    name="RMSE"
    )
PlotlyJS.plot([trace1]) # conicides with boxplot Figure 8


################## Generating backward values Table 1

size_historical_data = 16  # number of historical data to use: 8,14,15,16,30,60
t_k = 178 # (178 - Agust 15, 2020)
size_forecast = 30 # size_forecast
size_rolling_window = 0 
backward_case = Iran_backward(size_historical_data,t_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = backward_case[1] # approximation: 244.5951

#Table 1 (continuation)
#size_historical_data   RMSE_mean_estimated 
#8                      1269.3310
#12                     936.5389
#14                     405.2554
#15                     386.6915
#16                     286.3993
#17                     244.5951
#18                     3023.2656
#30                     704.9777
#60                     1470.2904

################## plot forecast forward and comparision MLP method

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

r = CSV.read("rateIran_paper.csv", DataFrame) # forecast paper until 15 sept - change Fecha
plotlyjs()

N = r.N[1] # N population size
I_scaled = N^(-1) * r.I0 #scaled infected data

d = CSV.read("dat_paper.csv", DataFrame)
aux2 = RSME(d.fore,I_scaled[(length(I_scaled) - size_forecast ):(end)]*N) # 407.60



p = Plots.plot(r.Fecha, r.I0,
    line=(2, :solid),
    label="Number infected per day", color="red",
    size=(675, 200), titleposition=:left,
    legend=:topright)


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
