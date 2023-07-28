# simulation of persistence examples (Figure 5)

include("./SDE_SIS.jl")
using Plots
using PlotlyJS
using Distributions

############################ INPUT ###########################################


################# sigma extremes
sigma_l = 0.07
sigma_u = 0.11
delta_sigma = sigma_u - sigma_l

################# beta extremes
beta_l = 0.25
beta_u = 0.30
delta_beta = beta_u - beta_l

################## TimeSubdivitions
Time_extrem = 300
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



################ Plot I_t evolution of 300 trajectories - Example Theorem 4
replications = 10000
gamma = 0.2 
mean = 0.0
I_0 = 0.01 # suppose that 30% of population are infected
delta = 0.001

Time_extrem = 900
subdivitions = 50
Time = range(0.001, Time_extrem, subdivitions)
h = 1 / subdivitions

################# Number deterministic basic reproduction R^D - example Theorem 4
R_D = (beta_l + beta_u) * 0.5 / gamma # 1.375
lim_It = (1 - 1 / R_D) # 0.2727

psi_ll = (1 / sigma_l^2) * (sqrt(beta_l^2 - 2 * sigma_l^2 * gamma) - (beta_l - sigma_l^2)) #Equation (31) Theorem 4
psi_lu = (1 / sigma_l^2) * (sqrt(beta_u^2 - 2 * sigma_l^2 * gamma) - (beta_u - sigma_l^2)) #Equation (32) Theorem 4


Trajectorie = zeros(replications, subdivitions)
M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
I = M[:, 2]

Trajectorie[1, :] = I
for i = 2:(replications-1)
    M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
    S = M[:, 1]
    I = M[:, 2]
    Trajectorie[i, :] = I

end
M = SDE_SIS(Time_extrem, subdivitions, h, mean, sigma, gamma, beta, I_0)
S = M[:, 1]
I = M[:, 2];
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
    label="Simulation I_t band", legend=:bottomright)
plot!(graf1, Time, psi_ll * ones(length(Time)), color="black", label="psi_l")#
plot!(graf1, Time, psi_lu * ones(length(Time)), color="black", label="psi_u")#

plot!(graf1, Time, lim_It * ones(length(Time)), color="red", label="lim It")#
