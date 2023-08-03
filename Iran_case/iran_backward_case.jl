
# -----------------------------
# -     BACKWARD
# -----------------------------

using Plots
using PlotlyJS
using DataFrames
using CSV
using Random, Distributions, LinearAlgebra
using Dates 
include("./Euler_Maruyama_SIS.jl") # only use SDE_SIS3


#Iran_forward(mj,j,size_forecast,size_rolling_window)

function Iran_backward(mj::Int64,j::Int64,size_forecast::Int64,size_rolling_window::Int64)
    # read historical data
    r = CSV.read("rateIran_paper.csv", DataFrame) # forecast paper until 15 sept - change Fecha
    plotlyjs()
    p = Plots.plot(r.Fecha, r.I0,
        line=(2, :solid),
        label="Number infected per day", color="red",
        size=(675, 200), titleposition=:left,
        legend=:topright)
    # scaling data
    N = r.N[1] # N population size
    I_scaled = N^(-1) * r.I0 #scaled data
    # define constants
    gamma = 1 / 14

    # time data subdivitions to use Euler Maruyama
    replications = 100000
    subdivitions = 50 # number of subdivitions into step values
    mean_noise = 0.0 
    alp = 0.05 # percentil to cut outiliers

    Trajectories = zeros(replications, size_forecast)

    initial_day_rolling_window = j - size_rolling_window
    rsme_data = zeros(length(I_scaled))
    for s in initial_day_rolling_window:(j)
        M = SDE_SIS5(subdivitions, size_forecast, mj, s, mean_noise, gamma, alp, I_scaled)
        Trajectories[1, :] = M[:, 2] # 1(suceptibles) 2(Infected)
        for i = 2:(replications-1)
            M = SDE_SIS5(subdivitions,  size_forecast, mj, s, mean_noise, gamma, alp, I_scaled)
            Trajectories[i, :] = M[:, 2]
        end
        M = SDE_SIS5(subdivitions, size_forecast, mj, s, mean_noise, gamma, alp, I_scaled)
        Trajectories[replications, :] = M[:, 2]
        Time = range(s, s + size_forecast - 1, size_forecast)
        mean_trajectories = zeros(length(Time))
        for i = 1:(length(Time))
            dat = Trajectories[:, i]
            mean_trajectories[i, 1] = quantile(dat, 0.5)#0.5
        end
        rsme_data[s] = RSME(mean_trajectories, I_scaled[s:s+size_forecast-1])
    end
    index_mj = mean(filter(!iszero, rsme_data))
    index_mj*N
    return [index_mj*N,Trajectories,p,rsme_data]  
end