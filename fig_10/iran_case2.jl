
# ----------------------------------
# -     Forecast SIS - data Iran 
# ----------------------------------

using Plots
using PlotlyJS
using DataFrames
using CSV
using Random, Distributions, LinearAlgebra
using Dates,Statistics
include("./Euler_Maruyama_SIS.jl") # use SDE_SIS3
include("./rmse.jl")

#Iran_forward(mj,j,size_forecast,size_rolling_window)

#size_forecast = 14 # testing
#size_day = 4 # testing
#mj= 60 # testing
#j = 75 # testing
#gamma = 1 / 14
#################### time data subdivitions to use Euler Maruyama
#replications = 1000 #testing
#subdivitions = 1000 # 50 # number of subdivitions into step values
#mean_noise = 0.0 # testing
#alp = 0.05 # percentil to cut outiliers

 ## testing
 #RSME(I_scaled[j+1:j + size_forecast],M[:,2])

function Iran_get_forecast_j(mj::Int64,j::Int64,size_forecast::Int64,size_day::Int64,replications::Int64,subdivitions::Int64,mean_noise::Float64,gamma::Float64)
    # read historical data
    r = CSV.read("rateIran_paper.csv", DataFrame) # forecast paper until 15 sept - change Fecha
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
    p = Plots.plot(r.Fecha, Float_I,
        line=(2, :solid),
        label="Number infected per day", color="red",
        size=(675, 200), titleposition=:left,
        legend=:topright,
        xticks=(r.Fecha, Dates.format.(r.Fecha, "yyyy-mm-dd")))
    # scaling data
    N = r.N[1] # N population size
    I_scaled = N^(-1) * Float_I #scaled data
    Trajectories = zeros(replications, size_forecast)

    #initial_day_rolling_window 
    rsme_data = zeros(length(I_scaled))
    #for s in initial_day_rolling_window:(j)
        #s = j#testing
        M = SDE_SIS1(subdivitions, size_forecast, size_day, mj, j, mean_noise, gamma, I_scaled)
        
        Trajectories[1, :] = M[:, 2] # 1(suceptibles) 2(Infected)
        for i = 2:(replications-1)
            M = SDE_SIS1(subdivitions, size_forecast, size_day, mj, j, mean_noise, gamma, I_scaled)
            Trajectories[i, :] = M[:, 2]
        end
        M = SDE_SIS1(subdivitions, size_forecast, size_day, mj, j, mean_noise, gamma, I_scaled)
        Trajectories[replications, :] = M[:, 2]

        Time = range(j, j + size_forecast - 1, size_forecast)
        mean_trajectories = zeros(length(Time))
        # estimator sigma and beta 

        for i = 1:(length(Time))
            dat = Trajectories[:, i]
            mean_trajectories[i, 1] = median(dat)#0.5
        end
        RMSE_error = RSME(mean_trajectories*N, I_scaled[j:j+size_forecast-1]*N)
    #end
    #index_mj = mean(filter(!iszero, rsme_data))
    #index_mj*N
    # retun index_mj*N,Trajectories,p,rsme_data
    return [Trajectories,mean_trajectories,RMSE_error,p]
end




