
using CSV
using Plots
using PlotlyJS
using DataFrames
using Dates
using FileIO
using ImageMagick
using Random, Distributions

include("./compute_AIC.jl")
include("./create_beta_estimator.jl")
include("./create_sigma_estimator.jl")
include("./rmse.jl")
include("./det_SIS.jl")

r = CSV.read("dat_gn.csv", DataFrame) # forecast paper until 15 sept - change Fecha

N = mean(r.Nt)  #r.N[1] # N population size
I_scaled = convert(Array{Float64},r.I_1) * N^(-1) 
    
    plotlyjs()
    # xticks=(r.Fecha, Dates.format.(r.Fecha, "yyyy-mm-dd"))
    p = Plots.plot(r.Year, r.I,
        line=(2, :solid),
        label="Number infected per day", color="red",
        size=(675, 200), titleposition=:left,
        legend=:topright)

#I_scaled = [0.01, 0.012, 0.013, 0.015, 0.016, 0.017, 0.020, 0.018, 0.019]  # ejemplo
h = 1.0
gamma = 1/55
mj = 2

# model 1: - all constants fixed, not estimation k=0 in compute_AIC
AIC_fixed_value = compute_AIC(I_scaled, h, gamma, mj)
println("fixed_AIC: ", AIC_fixed_value) # -228.94560610548947

# model 2: beta and sigma estimated fixed
AIC_deterministic_value = AIC_SIS_deterministic(I_scaled,mj, (length(I_scaled)-1))
println("det_AIC: ", AIC_deterministic_value) # -568.4538655710879

# model 3: - beta_t and sigma_t estimated
AIC_beta_sigma_value = compute_AIC(I_scaled, h, gamma, mj)
println("MLE_AIC: ", AIC_beta_sigma_value) # 87.05439389451053

# model 4: beta fixed and sigma_t estimated 
beta = 0.1 # for all t  
AIC_sigma_value = compute_AIC1(I_scaled, h, gamma,beta, mj)
println("MLE_AIC: ", AIC_sigma_value) # -72.94516538355998

# model 5: beta_t estimated and sigma_t fixed   
AIC_beta_value = compute_AIC2(I_scaled, h, gamma, mj)
println("MLE_AIC: ", AIC_beta_value) # -943.5276118768529


