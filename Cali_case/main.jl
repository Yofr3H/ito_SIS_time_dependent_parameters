using CSV
using Plots
using Plots
using PlotlyJS
using DataFrames

include("./cali_forward_case.jl")
include("./cali_backward_case.jl")

################## Generating forward values - Table 2
t_k = 563 # (563 - date: September 15, 2021 / 657 - date: December 18, 2021 )
size_historical_data = 7  # number of historical data to use: 5,6,7
size_forecast = 8 # size_forecast: 8,15,30
size_rolling_window =0 # 0: case no rolling window

forward_case = Cali_forward(size_historical_data,t_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = forward_case[1] # approximation: 12.1702 - Table 3
################### Table 2
#size_forecast  size_historical_data  RMSE_mean_estimated 	
#8		        5                     27.2238	             
#8		        6	                  23.2122	             
#8		        7	                  12.1702                 
#8		        8               	  25.5681              	
#15		        9                    344.7666	             
#15		        10	                  43.6809                
#15		        11	                  13.5403	             
#15             12                    23.6888	              
#30		        20	                  97.7435               
#30		        21	                  36.7797	              
#30		        22	                  81.4679

################## Generating backward values Table 3
t_k = 563 # (563 - date: September 15, 2021 )
size_historical_data = 7  # number of historical data to use: 5,6,7
size_forecast = 8 # size_forecast: 8,15,30
size_rolling_window =0 # 0: case no rolling window

backward_case = Cali_backward(size_historical_data,t_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = backward_case[1] # approximation: 15.7911 
#Table 4
#size forecast 	size_historical_data   	RMSE_mean_estimated 	
#8		        6	                    23.8016                
#8		        7                       11.9798	               
#8		        8	                    25.8656                
#15		        10	                    42.0432                
#15		        11	                    13.8436              
#15		        12                      24.5291 	               
#30		        20	                    97.7793                
#30		        21	                    90.6071	               
#30		        19                      95.4522	               

t_k = 657 # (657 - date: December 18, 2021 )
size_historical_data = 10  # number of historical data to use: 5,6,7
size_forecast = 15 # size_forecast: 8,15,30
size_rolling_window =0 # 0: case no rolling window

backward_case = Cali_backward(size_historical_data,t_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = backward_case[1]
forward_case = Cali_forward(size_historical_data,t_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = forward_case[1] # approximation: 12.9631 - Table 4
#Table 4
#size forecast	size_historical_data   	RMSE mean forward 	    RMSE mean backward
#8		        5                       861.5882	            873.7412 
#8		        6                       787.7186	            796.1732
#8		        7                       847.3904                853.0158
#15             8                       3070.9562               3074.5225
#15		        9                       1372.0081	            1367.6038
#15		        10                      1432.8608	            1427.7673	



Plot_original_data = forward_case[3] # Figure 9


################## plot forecast band

t_k = 657 # (563 - date: September 15, 2021 / 657 - date: December 18, 2021 )
size_historical_data = 6  # number of historical data to use: 5,6,7
size_forecast = 8 # size_forecast: 8,15,30
size_rolling_window =0 # 0: case no rolling window

forward_case = Cali_forward(size_historical_data,t_k,size_forecast,size_rolling_window)


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
data_boxplot_forward_case = Cali_forward(size_historical_data,t_k,size_forecast,size_rolling_window)

data_boxplot_forward = filter(!iszero, data_boxplot_forward_case[4]*N)
boxplot_forward8 = box(
    x=data_boxplot_forward,
    boxpoints="all",
    name="RMSE 8"
    )
    # run separately - preserve boxplot_forward# to plot
#= boxplot_forward30 = box(
    x=data_boxplot_forward,
    boxpoints="all",
    name="RMSE 30"
    ) =#
#= boxplot_forward15 = box(
    x=data_boxplot_forward,
    boxpoints="all",
    name="RMSE 15"
    ) =#

b = PlotlyJS.plot([boxplot_forward8, boxplot_forward15,boxplot_forward30]) # coincides with boxplot case a), Figure 10


################## Generating Forward values - Table 4

t_k = 657 # (657 - date: December 18, 2021 )
size_historical_data = 5  # number of historical data to use: 5,6,7
size_forecast = 8 # size_forecast: 8,15,30
size_rolling_window =0 # 0: case no rolling window

forward_case = Cali_forward(size_historical_data,t_k,size_forecast,size_rolling_window)
RMSE_mean_estimated = forward_case[1] # approximation: 15.7911 
#Table 4
#size forecast 	size_historical_data   	RMSE_mean_estimated 	
#8		        5                       916.6769	           
#8		        6                       890.2948	            
#8		        7                       895.8889                
#30		        19                      2196.5654	            	
#30		        20                      2199.3651              
#30		        21                      2195.6696	            s