
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

function gn_get_forecast(mj::Int64,j::Int64,size_forecast::Int64,size_day::Int64,replications::Int64,subdivitions::Int64,alp::Float64,mean_noise::Float64,gamma::Float64)
    # read historical data
    r = CSV.read("dat_gn.csv", DataFrame) # forecast paper until 15 sept - change Fecha
    ####################
    # scaling data
    N = mean(r.Nt)  #r.N[1] # N population size
    I_scaled = convert(Array{Float64},r.I_1) * N^(-1)  # r.I #scaled infected data
    

    # scaling data
    Trajectories = zeros(replications, size_forecast)

    #initial_day_rolling_window 
    rsme_data = zeros(length(I_scaled))
    #for s in initial_day_rolling_window:(j)
        #s = j#testing
        M = SDE_SISg(subdivitions, size_forecast, size_day, mj, j, mean_noise, gamma,alp, I_scaled)
        
        Trajectories[1, :] = M[:, 2] # 1(suceptibles) 2(Infected)
        for i = 2:(replications-1)
            M = SDE_SISg(subdivitions, size_forecast, size_day, mj, j, mean_noise, gamma,alp, I_scaled)
            Trajectories[i, :] = M[:, 2]
        end
        M = SDE_SISg(subdivitions, size_forecast, size_day, mj, j, mean_noise, gamma,alp, I_scaled)
        Trajectories[replications, :] = M[:, 2]

        Time = range(j, j + size_forecast - 1, size_forecast)
        mean_trajectories = zeros(length(Time))
        # estimator sigma and beta 

        for i = 1:(length(Time))
            dat = Trajectories[:, i]
            mean_trajectories[i, 1] = median(dat)#0.5
        end
        #RMSE_error = RSME(mean_trajectories*N, I_scaled[j:j+size_forecast-1]*N)
        RMSE_error = RSME(mean_trajectories, I_scaled[j:j+size_forecast-1])
    #end
    #index_mj = mean(filter(!iszero, rsme_data))
    #index_mj*N
    # retun index_mj*N,Trajectories,p,rsme_data
    return [Trajectories,mean_trajectories,RMSE_error]
end





size_forecast = 7 # testing
size_day = 2 # testing
mj= 20 # testing
j = 30 # testing
gamma = 1/55
#################### time data subdivitions to use Euler Maruyama
replications = 1000 #testing
subdivitions = 1000 # 50 # number of subdivitions into step values
mean_noise = 0.0 # testing
alp = 0.05 # percentil to cut outiliers


#M = SDE_SISg(subdivitions, size_forecast, size_day, mj, j, mean_noise, gamma, alp, I_scaled)
        

#model = gn_get_forecast(mj,j,size_forecast,size_day,replications,subdivitions,alp,mean_noise,gamma) # quantile procedure


