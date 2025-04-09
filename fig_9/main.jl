using CSV
using Plots
using PlotlyJS
using DataFrames
using Dates

include("./iran_case.jl")
include("./iran_case2.jl")

#------------------
# Initial parameters
#------------------

size_forecast = 15 # testing 15,30 
size_day = 2 # testing
mj= 168# testing 168
j = 178 # testing
gamma = 1 / 14
#################### time data subdivitions to use Euler Maruyama
replications = 10000 #testing
subdivitions = 1000 # 50 # number of subdivitions into step values
mean_noise = 0.0 # testing
alp = 0.05 

################## Generating Trajectories

model = Iran_get_forecast(mj,j,size_forecast,size_day,replications,subdivitions,alp,mean_noise,gamma) # quantile procedure
#model = Iran_get_forecast_j(mj,j,size_forecast,size_day,replications,subdivitions,alp,mean_noise,gamma) # final estimator

Trajectories = model[1]
mean_trajectories = model[2]
RMSE_mean_estimated = model[3]
Plot_original_data = model[4]

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

r = CSV.read("rateIran_paper.csv", DataFrame) # 
#r = CSV.read("rateCali.csv", DataFrame) # forecast paper until 15 sept - change Fecha
      
N = r.N[1] # N population size
I_scaled = N^(-1) * r.I0 #scaled infected data



plotlyjs()

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
     #xticks=(range(r.Fecha[mj],r.Fecha[end],8), "yy-mm-dd"))
     d = CSV.read("dat_paper.csv", DataFrame)
p = Plots.plot!(p,r.Fecha[178:208],d.fore,
    line=(2, :solid),
    label="MLP forecast", color="green",
    size=(675, 440), titleposition=:left) #case iran only

p = Plots.plot!(p,r.Fecha[j:(j+size_forecast-1)],median_trajectories[:, 1]*N,
    line=(2, :solid),
    label="Median forecast Itô-SIS", color="blue",
    size=(675, 440), titleposition=:left)
    #p = Plots.plot!(p,r.Fecha[j:(j+size_forecast-1)],mean_trajectories[:, 1]*N,
    #line=(2, :solid),
    #label="Mean forecast Itô-SIS", color="orange",
    #size=(675, 200), titleposition=:left,
    #legend=:bottomright)
p = Plots.plot(p,r.Fecha[j:(j+size_forecast-1)], extrem_trajectories[:, 1]*N,
    fillrange=extrem_trajectories[:, 2]*N, fillalpha=0.15, c=1,
    label="forecast Itô-SIS band")
p
Plots.savefig("p.svg")

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
r = CSV.read("rateIran_paper.csv", DataFrame) #           
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
        #Plot_original_data = model[4]

        #Time = range(j, j + size_forecast - 1, size_forecast)
        #extrem_trajectories = zeros(length(Time), 2)
        #median_trajectories = zeros(length(Time))
        #for i = 1:(length(Time))
            dat = Trajectories[:, 2]
            mean_forecast[j, 1] = mean(dat) #quantile(dat, 0.5)#0.95
            extrem_forecast[j, 1] = minimum(dat)
            extrem_forecast[j, 2] = maximum(dat)
        #end
end

r = CSV.read("rateIran_paper.csv", DataFrame) # forecast paper until 15 sept - change Fecha
#r = CSV.read("rateCali.csv", DataFrame) # forecast paper until 15 sept - change Fecha
 
#################### eliminate null Infected values
Float_I = Float64.(r.I0) 
for k in 1:length(r.I0)
    if Float_I[k] == 0
        Float_I[k] = 0.5 * (Float_I[k-1] + Float_I[k+1] )
    end
 end
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

#################### run only  iran case
#size_forecast = 30
d = CSV.read("dat_paper.csv", DataFrame)
#aux2 = RSME(d.fore[1:(size_forecast+1)],I_scaled[(length(I_scaled) - size_forecast ):(end)]*N) # 407.60

g = Plots.plot!(g,r.Fecha[178:208],d.fore,
    line=(2, :solid),
    label="MLP forecast", color="green",
    size=(1495, 640), titleposition=:left) #case iran only
#################################

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
Plots.savefig("g_rolling.svg")
 g2=g   
Plots.savefig("g2.png") # fig 6



















 # approximation: 386.13356

#Table 1 - Iran case
#size_historical_data   RMSE_mean_estimated 
#6                      1041.7964
#7                      1158.5953
#8                      1231.74026
#9                      1012.0450
#10                     860.65037
#11                     595.1842
#12                     386.1335
#13                     1152.7663
#29                     570.2873
#30                     713.4612
#31                     454.8572
#32                     1518.6082
#60                     1470.0902



