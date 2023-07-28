#
# From the estimation of parameters β γ and σ (forward approximation), 
# a discretization of Itô-SIS.
# Simulate empirical distribution of β γ and σ with hommogenized data
# free of missing data or outiliers.
#

# -----------------------------
# -     FORWARD CASE
# -----------------------------

using Random, Distributions, StatsBase

include("./create_forward_beta_estimator.jl")
include("./create_forward_sigma_estimator.jl")
include("./create_backward_beta_estimator.jl")
include("./create_backward_sigma_estimator.jl")
include("./rmse.jl")


function SDE_SIS3(subdivitions::Int64, size_forecast::Int64, mj::Int64, j::Int64, mean_noise::Float64, gamma::Float64, alp::Float64, I_scaled::Any)
    T = subdivitions
    Time = size_forecast
    h = 1 / T
    desvest = zeros(length(I_scaled))
    k = j - mj + 1 
    for i in k:j
        desvest[i] = create_sigma_estimator(I_scaled[1:(i)], mj, i)
    end

    beta = zeros(length(I_scaled))
    for i in k:j
        beta[i] = create_beta_estimator(I_scaled[1:i], gamma, mj, i)
    end
    S1 = zeros(T + 1, 2)#auxiliar intermediate steps: 1 (Suceptible) 2(Infected)
    #create_vectors of returns STEP 1
    S = zeros(Time)
    I = zeros(Time)
    I0 = I_scaled[j]
    #Initialization of I_0 S_0 
    S1[1, 2] = I0 #I_0
    S1[1, 1] = 1.0 - I0 #S_0
    S[1] = S1[1, 1]
    I[1] = S1[1, 2]
    # SIMULATION STEP 2
    for s in 2:(T+1)
        rand_clean = rand(Uniform(alp, 1 - alp))
        desv_sim = quantile(desvest[k:j], rand_clean)
        z = rand(Normal(mean_noise, 1.0))
        #evaluate Infected next  step - rate infection with additive noise
        bt = quantile(beta[k:j], rand_clean)
        S1[s, 2] = S1[s-1, 2] + (bt + desv_sim*z*sqrt(h)) * (1.0 - S1[s-1, 2]) * S1[s-1, 2] * h - gamma * S1[s-1, 2] * h#Infected
        S1[s, 1] = 1.0 - S1[s, 2]
    end
    #saving evolution of second day (day is divided in T parts - T+1 extrems)
    S[2] = S1[T+1, 1]
    I[2] = S1[T+1, 2]

    # Use forecasting to improve the empirical distribution
    new_I_scaled = zeros(length(I_scaled)) #historical data and forecasting
    new_I_scaled[1:j] = I_scaled[1:j]
    new_I_scaled[j+1] = I[2]#obtain new estimators of sigma and beta
    #obtain new estimators of sigma and beta
    desvest[j+1] = create_sigma_estimator(new_I_scaled[1:(j+1)], mj, (j + 1))
    beta[j+1] = create_beta_estimator(new_I_scaled[1:(j+1)], gamma, mj, (j + 1))

    #updating to recalculate next iteration
    S1 = zeros(T + 1, 2)
    S1[1, 1] = S[2]
    S1[1, 2] = I[2]
  
    for i in 3:Time
        for w in 2:(T+1)
            rand_clean = rand(Uniform(alp, 1 - alp))
            desv_sim = quantile(desvest[k:(j+(i-2))], rand_clean)
            bt = quantile(beta[k:(j+(i-2))], rand_clean)
            z = rand(Normal(mean_noise, 1.0))
            #evaluate Infected next  step - rate infection with additive noise
            S1[w, 2] = S1[w-1, 2] + (bt + desv_sim*z*sqrt(h)) * (1.0 - S1[w-1, 2]) * S1[w-1, 2] * h - gamma * S1[w-1, 2] * h#Infected
            #evaluate Infected next  step - rate infection with proportional noise
            S1[w, 1] = 1.0 - S1[w, 2]
        end
        #saving
        S[i] = S1[T+1, 1]
        I[i] = S1[T+1, 2]

        
        new_I_scaled[j+(i-1)] = I[i]#obtain new estimators of sigma and beta
        #obtain new estimators of sigma and beta
        desvest[j+(i-1)] = create_sigma_estimator(new_I_scaled[1:(j+(i-1))], mj, (j + (i-1) ))
        beta[j+(i-1)] = create_beta_estimator(new_I_scaled[1:(j+(i-1))], gamma, mj, (j + (i-1)))
     
        #updating values
        S1 = zeros(T + 1, 3)
        S1[1, 1] = S[i]
        S1[1, 2] = I[i]
    end
    return [S I]
end


function SDE_SIS5(subdivitions::Int64, size_forecast::Int64, mj::Int64, j::Int64, mean_noise::Float64, gamma::Float64, alp::Float64, I_scaled::Any)
    # use create_backward_beta_estimator and create_backward_sigma_estimator
     T = subdivitions
     Time = size_forecast
     h = 1 / T
     desvest = zeros(length(I_scaled))
     k = j - mj + 1 #
     for i in k:j
         desvest[i] = create_backward_sigma_estimator(I_scaled[1:(i)], mj, i)
     end
 
     beta = zeros(length(I_scaled))
 
     for i in k:j
         beta[i] = create_backward_beta_estimator(I_scaled[1:i], gamma, mj, i)
     end
 
     S1 = zeros(T + 1, 2)#auxiliar intermediate steps: 1 (Suceptible) 2(Infected)
     #create_vectors of returns STEP 1
     S = zeros(Time)
     I = zeros(Time)
     C = zeros(Time)
     I0 = I_scaled[j]
     #Initialization of I_0 S_0 
     S1[1, 2] = I0 #I_0
     S1[1, 1] = 1.0 - I0 #S_0
     S[1] = S1[1, 1]
     I[1] = S1[1, 2]
     # SIMULATION STEP 2
     for s in 2:(T+1)
         rand_clean = rand(Uniform(alp, 1 - alp))
         desv_sim = quantile(desvest[k:j], rand_clean)
         z = rand(Normal(mean_noise, 1.0))
         #evaluate Infected next  step - rate infection with additive noise
         rand_clean = rand(Uniform(alp, 1 - alp))
         bt = quantile(beta[k:j], rand_clean)
         S1[s, 2] = S1[s-1, 2] + (bt + desv_sim*z*sqrt(h)) * (1.0 - S1[s-1, 2]) * S1[s-1, 2] * h - gamma * S1[s-1, 2] * h#Infected
         S1[s, 1] = 1.0 - S1[s, 2]
     end
     #saving evolution of second day (day is divided in T parts - T+1 extrems)
     S[2] = S1[T+1, 1]
     I[2] = S1[T+1, 2]
 
     # idea: use forecasting to improve the empirical distribution
     new_I_scaled = zeros(length(I_scaled)) #historical data and forecasting
     new_I_scaled[1:j] = I_scaled[1:j]
     new_I_scaled[j+1] = I[2]#obtain new estimators of sigma and beta
     #obtain new estimators of sigma and beta
     desvest[j+1] = create_backward_sigma_estimator(new_I_scaled[1:(j+1)], mj, (j + 1))
     beta[j+1] = create_backward_beta_estimator(new_I_scaled[1:(j+1)], gamma, mj, (j + 1))
 
     #updating to recalculate next iteration
     S1 = zeros(T + 1, 2)
     S1[1, 1] = S[2]
     S1[1, 2] = I[2]
 
   
     for i in 3:Time
         for w in 2:(T+1)
             rand_clean = rand(Uniform(alp, 1 - alp))
             desv_sim = quantile(desvest[k:(j+(i-2))], rand_clean)
             bt = quantile(beta[k:(j+(i-2))], rand_clean)
             z = rand(Normal(mean_noise, 1.0))
             #evaluate Infected next  step - rate infection with additive noise
             S1[w, 2] = S1[w-1, 2] + (bt + desv_sim*z*sqrt(h)) * (1.0 - S1[w-1, 2]) * S1[w-1, 2] * h - gamma * S1[w-1, 2] * h#Infected
             #evaluate Infected next  step - rate infection with proportional noise
             S1[w, 1] = 1.0 - S1[w, 2]
         end
         #saving
         S[i] = S1[T+1, 1]
         I[i] = S1[T+1, 2]
 
         
         new_I_scaled[j+i] = I[i]#obtain new estimators of sigma and beta
         #obtain new estimators of sigma and beta
         desvest[j+(i-1)] = create_sigma_estimator(new_I_scaled[1:(j+(i-1))], mj, (j + (i-1) ))
         beta[j+(i-1)] = create_beta_estimator(new_I_scaled[1:(j+(i-1))], gamma, mj, (j + (i-1)))
 
         #updating values
         S1 = zeros(T + 1, 3)
         S1[1, 1] = S[i]
         S1[1, 2] = I[i]
     end
     return [S I]
 end
 