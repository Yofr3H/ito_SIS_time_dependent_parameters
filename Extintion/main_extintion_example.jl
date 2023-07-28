# simulation of extintion example (Figures 1,2,3 and 4)

include("./SDE_SIS.jl")
using Plots
using PlotlyJS
using Distributions

############################ INPUT ###########################################


################# sigma extremes
sigma_l = 0.65
sigma_u = 0.75
delta_sigma = sigma_u - sigma_l

################# beta extremes
beta_l = 0.45
beta_u = 0.50
delta_beta = beta_u - beta_l

################## TimeSubdivitions
Time_extrem = 400
subdivitions = 50

Time = range(0.001, Time_extrem, subdivitions)
h = 1 / subdivitions

################## beta function
y = range(0.0, 1 * Time_extrem, subdivitions)
beta = (delta_beta) * 0.25 * sin.(pi * 2.0 *  0.3 * y) .+ delta_beta * 0.25 * cos.( pi * 2.0 * (delta_beta) * 0.3 * y) .+ beta_l .+ delta_beta * 0.5
plotlyjs()
Plots.plot(y, beta, color="red", label="β_t")

################# sigma function
x = range(0.0, 1 * Time_extrem, subdivitions)
sigma = (delta_sigma) * 0.25 * sin.(pi * 2.0 * 0.3275 * x) .+ delta_sigma * 0.25 * cos.(pi * 2.0 * (delta_sigma) * 0.03275 * x) .+ sigma_l .+ delta_sigma * 0.5
plotlyjs()
Plots.plot(x, sigma, color="blue", label="σ_t")




################ Plot I_t evolution of 500 trajectories - Example Theorem 2
replications = 500
gamma = 0.7 
mean = 0.0
I_0 = 0.3 # suppose that 30% of population are infected
delta = 0.001

Time_extrem = 40
subdivitions = 50
Time = range(0.001, Time_extrem, subdivitions)
h = 1 / subdivitions

Trajectorie = zeros(replications, subdivitions)
M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
I = M[:, 2]
Trajectorie[1, :] = I
for i = 2:(replications-1)
    M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
    S = M[:, 1]
    I = M[:, 2]
    Trajectorie[i, :] = I#I
end
M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
S = M[:, 1]
I = M[:, 2]
Trajectorie[replications, :] = I#I
extrem_trajectorie = zeros(length(Time), 2)
mean_trajectorie = zeros(length(Time))
for i = 1:(length(Time))
    dat = Trajectorie[:, i]
    mean_trajectorie[i, 1] = quantile(dat, 0.5)
    extrem_trajectorie[i, 1] = minimum(dat)
    extrem_trajectorie[i, 2] = maximum(dat)
end
graf1 = Plots.plot(Time, extrem_trajectorie[:, 1],
    fillrange=extrem_trajectorie[:, 2], fillalpha=0.15, c=1,
    label="Simulation I_t band", legend=:topleft)

################ Plot log(I_t)/t evolution of 500 trajectories - Example Theorem 2
replications = 500
gamma = 0.7 
mean = 0.0
I_0 = 0.3 # suppose that 30% of population are infected
delta = 0.001

Time_extrem = 40
subdivitions = 50
Time = range(0.001, Time_extrem, subdivitions)
h = 1 / subdivitions

Trajectorie = zeros(replications, subdivitions)
M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
I = M[:, 2]
I = log.(I)./Time#prove extintion theorem
Trajectorie[1, :] = I
for i = 2:(replications-1)
    M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
    S = M[:, 1]
    I = M[:, 2]
    I = log.(I)./Time #prove extintion theorem
    Trajectorie[i, :] = I

end
M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
S = M[:, 1]
I = M[:, 2];
I = log.(I) ./ Time#prove extintion theorem
Trajectorie[replications, :] = I#I
extrem_trajectorie = zeros(length(Time), 2)
mean_trajectorie = zeros(length(Time))
for i = 1:(length(Time))
    dat = Trajectorie[:, i]
    mean_trajectorie[i, 1] = quantile(dat, 0.5)#0.95
    extrem_trajectorie[i, 1] = minimum(dat)
    extrem_trajectorie[i, 2] = maximum(dat)
end
graf2 = Plots.plot(Time, extrem_trajectorie[:, 1],
    fillrange=extrem_trajectorie[:, 2], fillalpha=0.15, c=1,
    label="Simulation log(I_t)/t band", legend=:topright)

################## Simulation example Theorem 3

################# sigma extremes
sigma_l = 0.72
sigma_u = 0.78
delta_sigma = sigma_u - sigma_l

################# sigma function
x = range(0.0, 1 * Time_extrem, subdivitions)
sigma = (delta_sigma) * 0.25 * sin.(pi * 2.0 * 0.3275 * x) .+ delta_sigma * 0.25 * cos.(pi * 2.0 * (delta_sigma) * 0.03275 * x) .+ sigma_l .+ delta_sigma * 0.5
plotlyjs()
Plots.plot(x, sigma, color="blue", label="σ_t")


################ Plot I_t evolution of 500 trajectories - Example Theorem 3
replications = 500
gamma = 0.7 
mean = 0.0
I_0 = 0.3 # suppose that 30% of population are infected
delta = 0.001

Time_extrem = 40
subdivitions = 50
Time = range(0.001, Time_extrem, subdivitions)
h = 1 / subdivitions

Trajectorie = zeros(replications, subdivitions)
M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
I = M[:, 2]
Trajectorie[1, :] = I
for i = 2:(replications-1)
    M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
    S = M[:, 1]
    I = M[:, 2]
    #   I = log.(I)./Time#prove extintion theorem
    Trajectorie[i, :] = I#I
    #plot!(graf, Time, I, color="green", title="", label="")#I

end
M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
S = M[:, 1]
I = M[:, 2];
#I = log.(I) ./ Time#prove extintion theorem
Trajectorie[replications, :] = I#I
extrem_trajectorie = zeros(length(Time), 2)
mean_trajectorie = zeros(length(Time))
for i = 1:(length(Time))
    dat = Trajectorie[:, i]
    mean_trajectorie[i, 1] = quantile(dat, 0.5)#0.95
    extrem_trajectorie[i, 1] = minimum(dat)
    extrem_trajectorie[i, 2] = maximum(dat)
end
graf1 = Plots.plot(Time, extrem_trajectorie[:, 1],
    fillrange=extrem_trajectorie[:, 2], fillalpha=0.15, c=1,
    label="Simulation I_t band", legend=:topleft)

################ Plot log(I_t)/t evolution of 500 trajectories - Example Theorem 3
replications = 500
gamma = 0.7 
mean = 0.0
I_0 = 0.3 # suppose that 30% of population are infected
delta = 0.001

Time_extrem = 40
subdivitions = 50
Time = range(0.001, Time_extrem, subdivitions)
h = 1 / subdivitions

Trajectorie = zeros(replications, subdivitions)
M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
I = M[:, 2]
I = log.(I)./Time#prove extintion theorem
Trajectorie[1, :] = I
for i = 2:(replications-1)
    M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
    S = M[:, 1]
    I = M[:, 2]
    I = log.(I)./Time #prove extintion theorem
    Trajectorie[i, :] = I

end
M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
S = M[:, 1]
I = M[:, 2];
I = log.(I) ./ Time#prove extintion theorem
Trajectorie[replications, :] = I#I
extrem_trajectorie = zeros(length(Time), 2)
mean_trajectorie = zeros(length(Time))
for i = 1:(length(Time))
    dat = Trajectorie[:, i]
    mean_trajectorie[i, 1] = quantile(dat, 0.5)#0.95
    extrem_trajectorie[i, 1] = minimum(dat)
    extrem_trajectorie[i, 2] = maximum(dat)
end
graf2 = Plots.plot(Time, extrem_trajectorie[:, 1],
    fillrange=extrem_trajectorie[:, 2], fillalpha=0.15, c=1,
    label="Simulation log(I_t)/t band", legend=:topright)
