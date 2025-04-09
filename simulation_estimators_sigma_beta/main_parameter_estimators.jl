# simulation to show performance Unbiased estimator \sigma(t) and MLE(\beta)(t)


using Plots
using PlotlyJS
using Distributions
using FileIO
using ImageMagick
include("./SDE_SIS.jl")

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
Time_extrem = 200
divitions = 300
subdivitions = 50
replications = 10000 #replications

Time = range(0.00, Time_extrem, divitions)


################## beta function
y = range(0.0, 1 * Time_extrem, subdivitions * divitions)
beta = (delta_beta) * 0.25 * sin.(pi * 2.0 *  0.3 * y) .+ delta_beta * 0.25 * cos.( pi * 2.0 * (delta_beta) * 0.3 * y) .+ beta_l .+ delta_beta * 0.5
plotlyjs()
beta_figure_consistency = Plots.plot(y, beta, color="red", label="β_t")
Plots.savefig("beta_figure_consistency.png")
img33 = load("beta_figure_consistency.png")
save("beta_figure_consistency.jpg",img33)

################# sigma function
x = range(0.0, 1 * Time_extrem, subdivitions * divitions)
sigma = (delta_sigma) * 0.25 * sin.(pi * 2.0 * 0.3275 * x) .+ delta_sigma * 0.25 * cos.(pi * 2.0 * (delta_sigma) * 0.03275 * x) .+ sigma_l .+ delta_sigma * 0.5
plotlyjs()
sigma_figure_consistency = Plots.plot(x, sigma, color="blue", label="σ_t")
Plots.savefig("sigma_figure_consistency.png")
img13 = load("sigma_figure_consistency.png")
save("sigma_figure_consistency.jpg",img13)
#################### strong consistency for beta and sigma #####################
gamma = 0.7 
mean1 = 0.0
I_0 = 0.3 # suppose that 30% of population are infected


beta_trajectorie = zeros(replications, length(Time))
sigma_trajectorie = zeros(replications, length(Time))
for i = 2:(replications-1)
    M = SDE_SIS(Time_extrem, divitions, subdivitions, mean1, sigma, gamma, beta, I_0)
    beta_trajectorie[i, :] = M[:,3]
    sigma_trajectorie[i, :] = M[:,4]
end

mean_MLE_beta_estimator = zeros(length(Time))
mean_sigma_estimator = zeros(length(Time))
for i = 1:(length(Time))
    aux1 = beta_trajectorie[:, i]
    aux2 = sigma_trajectorie[:, i]
    mean_MLE_beta_estimator[i, 1] = mean(aux1)
    mean_sigma_estimator[i, 1] = mean(aux2)
end

plotlyjs()

beta_estimator_figure = Plots.plot(Time, mean_MLE_beta_estimator, color="red", label="MLE_beta")
Plots.savefig("beta_estimator_figure.png")
img1313 = load("beta_estimator_figure.png")
save("beta_estimator_figure.jpg",img1313)
plotlyjs()
beta_estimator = Plots.plot(Time, mean_sigma_estimator, color="blue", label="hat(sigma)")
Plots.savefig("sigma_estimator_figure.png")
img313 = load("sigma_estimator_figure.png")
save("sigma_estimator_figure.jpg",img313)


function MAPE(y_true::Vector{Float64}, y_pred::Vector{Float64})
    if length(y_true) != length(y_pred)
        error("Vectors must have the same length")
    end
    return 100 * mean(abs.((y_true .- y_pred) ./ y_true))
end

Percentual_error_beta = MAPE(beta[1:50:end],mean_MLE_beta_estimator)

Percentual_error_sigma = MAPE(sigma[1:50:end],mean_sigma_estimator)

