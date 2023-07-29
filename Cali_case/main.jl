using CSV
using Plots
using Plots
using PlotlyJS
using DataFrames

include("./cali_forward_case.jl")
include("./cali_backward_case.jl")

################## Generating forward values - Table 2
I_k = 563 # (563 - date: September 15, 2021 / 657 - date: December 18, 2021 )
size_historical_data = 6  # number of historical data to use: 5,6,7
size_forecast = 8 # size_forecast: 8,15,30
size_rolling_window =0 # 0: case no rolling window

forward_case = Cali_forward(size_historical_data,I_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = forward_case[1] # approximation: 12.9631 - Table 2
################### Table 2
#size_forecast  size_historical_data  RMSE_mean_estimated 	% population
#8		        5                     12.9644	              0.0005136
#8		        6	                  10.5747	              0.0004189
#8		        7	                  11.9399                 0.0004730
#8		        8               	  12.0344                 0.0004768	
#15		        9                     17.6087	              0.0006976
#15		        10	                  14.3588                 0.0005689
#15		        11	                  18.8530	              0.0007469
#30		        19                    42.2818	              0.0016752
#30		        20	                  32.9803                 0.0013067
#30		        21	                  40.9528	              0.0016225

################## Generating backward values Table 3
I_k = 563 # (563 - date: September 15, 2021 )
size_historical_data = 5  # number of historical data to use: 5,6,7
size_forecast = 8 # size_forecast: 8,15,30
size_rolling_window =0 # 0: case no rolling window

backward_case = Cali_backward(size_historical_data,I_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = backward_case[1] # approximation: 15.7911 
#Table 3
#size forecast 	size_historical_data   	RMSE_mean_estimated 	% population
#8		        6	                    11.8185                 0.0004682
#8		        5                       15.7911	                0.0006256
#8		        7	                    13.1992                 0.0005229	
#15		        10	                    17.6456                 0.0006991
#15		        11	                    22.2810	                0.0008827
#15		        9                       15.6998	                0.0006220
#30		        20	                    35.1036                 0.0013908
#30		        21	                    42.7376	                0.0016933
#30		        19                      44.5962	                0.0017669

I_k = 657 # (657 - date: December 18, 2021 )
size_historical_data = 5  # number of historical data to use: 5,6,7
size_forecast = 8 # size_forecast: 8,15,30
size_rolling_window =0 # 0: case no rolling window

forward_case = Cali_forward(size_historical_data,I_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = forward_case[1] # approximation: 12.9631 - Table 4
#Table 4
#size forecast	size_historical_data   	RMSE_mean_estimated 	% population
#8		        5                       916.6769	            0.0363194
#8		        6                       890.2948	            0.0352742
#8		        7                       895.8889                0.0354958
#30		        19                      2196.5654	            0.0870297	
#30		        20                      2199.3651               0.0871406
#30		        21                      2195.6696	            0.0869942


Plot_original_data = forward_case[3] # Figure 9


################## plot forecast band

I_k = 657 # (563 - date: September 15, 2021 / 657 - date: December 18, 2021 )
size_historical_data = 6  # number of historical data to use: 5,6,7
size_forecast = 8 # size_forecast: 8,15,30
size_rolling_window =0 # 0: case no rolling window

forward_case = Cali_forward(size_historical_data,I_k,size_forecast,size_rolling_window)


Trajectories = forward_case[2]
j = I_k
Time = range(j, j + size_forecast - 1, size_forecast)
extrem_trajectories = zeros(length(Time), 2)
median_trajectories = zeros(length(Time))
for i = 1:(length(Time))
    dat = Trajectories[:, i]
    median_trajectories[i, 1] = quantile(dat, 0.5)#0.95
    extrem_trajectories[i, 1] = minimum(dat)
    extrem_trajectories[i, 2] = maximum(dat)
end

r = CSV.read("rateCali1.csv", DataFrame) 
plotlyjs()
p = Plots.plot(r.Fecha, r.I0,
    line=(2, :solid),
    label="Number infected per day", color="red",
    size=(675, 200), titleposition=:left,
    legend=:topright)

N = r.N[1] # N population size
I_scaled = N^(-1) * r.I0 #scaled infected data
p = Plots.plot!(p,r.Fecha[j:(j+size_forecast-1)],median_trajectories[:, 1]*N,
    line=(2, :solid),
    label="Median forecast Itô-SIS", color="blue",
    size=(675, 200), titleposition=:left,
    legend=:bottomright)
p = Plots.plot(p,r.Fecha[j:(j+size_forecast-1)], extrem_trajectories[:, 1]*N,
    fillrange=extrem_trajectories[:, 2]*N, fillalpha=0.15, c=1,
    label="forecast Itô-SIS band", legend=:topleft)
p # To see the band, Zoom near of 657 - date: December 18, 2021

################ Boxplot of RMSE median, using the rolling window - Figure 10 

size_rolling_window = 8 # case b) 15 and case c) 30
data_boxplot_forward_case = Cali_forward(size_historical_data,I_k,size_forecast,size_rolling_window)

data_boxplot_forward = filter(!iszero, data_boxplot_forward_case[4]*N)
boxplot_forward = box(
    x=data_boxplot_forward,
    boxpoints="all",
    name="RMSE"
    )
PlotlyJS.plot([boxplot_forward]) # coincides with boxplot case a), Figure 10


################## Generating Forward values - Table 4

I_k = 657 # (657 - date: December 18, 2021 )
size_historical_data = 5  # number of historical data to use: 5,6,7
size_forecast = 8 # size_forecast: 8,15,30
size_rolling_window =0 # 0: case no rolling window

forward_case = Cali_forward(size_historical_data,I_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = forward_case[1] # approximation: 15.7911 
#Table 4
#size forecast 	size_historical_data   	RMSE_mean_estimated 	% population
#8		        5                       916.6769	            0.0363194
#8		        6                       890.2948	            0.0352742
#8		        7                       895.8889                0.0354958
#30		        19                      2196.5654	            0.0870297	
#30		        20                      2199.3651               0.0871406
#30		        21                      2195.6696	            0.0869942