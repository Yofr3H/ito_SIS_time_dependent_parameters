#
# From the estimation of parameters β γ and σ (forward approximation), 
# a discretization of Itô-SIS.
# Simulate empirical distribution of β γ and σ with hommogenized data
# free of missing data or outiliers.
#

# -----------------------------
# -     FORWARD CASE
# -----------------------------

using Random, Distributions
#StatsBase

include("./create_beta_estimator.jl")
include("./create_sigma_estimator.jl")
include("./rmse.jl")


function SDE_SIS3(subdivitions::Int64, size_forecast::Int64,size_day::Int64, mj::Int64, j::Int64, mean_noise::Float64, gamma::Float64, alp::Float64, I_scaled::Any)
    T = subdivitions # = 1000 
    Time_size = size_forecast # = 14
    #size_day = 5
    desvest = zeros(length(I_scaled))
    Time_interval = range(mj, j, j-mj) 
    k = mj + 1
    h_data = Time_interval[2] - Time_interval[1] 
    for i in k:j
        desvest[i] = create_sigma_estimator(I_scaled,h_data, mj, i)
    end

    beta = zeros(length(I_scaled))
    for i in k:j
        beta[i] = create_beta_estimator(I_scaled,h_data, gamma, mj, i)
    end
    S1 = zeros(T + 1, 2)#auxiliar intermediate steps: 1 (Suceptible) 2(Infected)
    #create_vectors of returns STEP 1
    S = zeros(Time_size)
    I = zeros(Time_size)
    I0 = I_scaled[j]
    #Initialization of I_0 S_0 
    S1[1, 2] = I0 #I_0
    S1[1, 1] = 1.0 - I0 #S_0
    S[1] = S1[1, 1]
    I[1] = S1[1, 2]
    h =  T^(-1)
    # SIMULATION STEP 2

    for s in 2:(T+1)
        rand_clean = rand(Uniform(alp, 1 - alp)) # desvest[j] #quantile(desvest[j-size_forecast:j], rand_clean)#
        desv_sim =  quantile(desvest[k:j], rand_clean)
        z = rand(Normal(mean_noise, 1.0))
        #evaluate Infected next  step - rate infection with additive noise# beta[j] #quantile(beta[j - size_forecast:j], rand_clean)
        bt = quantile(beta[k:j], rand_clean)#
        S1[s, 2] = S1[s-1, 2] + (bt * h + desv_sim * z * sqrt(h)) * (1.0 - S1[s-1, 2]) * S1[s-1, 2]  - gamma * S1[s-1, 2] * h#Infected
        S1[s, 1] = 1.0 - S1[s, 2]
    end
    #saving evolution of second day (day was divided in T parts - T+1 extrems)
    I[2] = S1[(size_day + 1),2] # the forecast by day, the value 10
    S[2] = 1 - I[2]

    # Use forecasting to improve the empirical distribution
    new_I_scaled = zeros(length(I_scaled)) #historical data and forecasting
    new_I_scaled[1:j] = I_scaled[1:j]
    new_I_scaled[j+1] = I[2]#obtain new estimators of sigma and beta
    #obtain new estimators of sigma and beta
    desvest[j+1] = create_sigma_estimator(new_I_scaled[1:(j+1)],h_data, mj, (j + 1))
    beta[j+1] = create_beta_estimator(new_I_scaled[1:(j+1)],h_data, gamma, mj, (j + 1))

    #updating to recalculate next iteration
    S1 = zeros(T + 1, 2)
    S1[1, 1] = S[2]
    S1[1, 2] = I[2]
  
    for i in 3:Time_size
        for w in 2:(T+1) # quantile(beta[j - size_forecast:(j+(i-2))], rand_clean) #
            rand_clean = rand(Uniform(alp, 1 - alp))#   beta[j + i - 2]#  desvest[j + i - 2]# quantile(desvest[j - size_forecast:(j+(i-2))], rand_clean) #
            desv_sim = quantile(desvest[k:(j+(i-2))], rand_clean)# 
            bt = quantile(beta[k:(j+(i-2))], rand_clean)#
            z = rand(Normal(mean_noise, 1.0))
            #evaluate Infected next  step - rate infection with additive noise
            S1[w, 2] = S1[w-1, 2] + (bt* h  + desv_sim*z*sqrt(h)) * (1.0 - S1[w-1, 2]) * S1[w-1, 2] - gamma * S1[w-1, 2] * h#Infected
            #evaluate Infected next  step - rate infection with proportional noise
            S1[w, 1] = 1.0 - S1[w, 2]
        end
        #saving
        I[i] = S1[size_day + 1, 2] # save the first day, with position 10
        S[i] = 1 - I[i]

        
        new_I_scaled[j+(i-1)] = I[i]#obtain new estimators of sigma and beta
        #obtain new estimators of sigma and beta
        desvest[j+(i-1)] = create_sigma_estimator(new_I_scaled[1:(j+(i-1))],h, k, (j + (i-1) ))
        beta[j+(i-1)] = create_beta_estimator(new_I_scaled[1:(j+(i-1))],h, gamma, k, (j + (i-1)))
     
        #updating values
        S1 = zeros(T + 1, 3)
        S1[1, 1] = S[i]
        S1[1, 2] = I[i]
    end
    return [S I ]
end

# not use quatile only the final beta and sigma estimators
function SDE_SIS1(subdivitions::Int64, size_forecast::Int64,size_day::Int64, mj::Int64, j::Int64, mean_noise::Float64, gamma::Float64, alp::Float64, I_scaled::Any)
    T = subdivitions # = 1000 
    Time_size = size_forecast # = 14
    #size_day = 5
    desvest = zeros(length(I_scaled))
    Time_interval = range(mj, j, j-mj) 
    k = mj + 1
    h_data = Time_interval[2] - Time_interval[1] 
    for i in k:j
        desvest[i] = create_sigma_estimator(I_scaled,h_data, mj, i)
    end

    beta = zeros(length(I_scaled))
    for i in k:j
        beta[i] = create_beta_estimator(I_scaled,h_data, gamma, mj, i)
    end
    S1 = zeros(T + 1, 2)#auxiliar intermediate steps: 1 (Suceptible) 2(Infected)
    #create_vectors of returns STEP 1
    S = zeros(Time_size)
    I = zeros(Time_size)
    I0 = I_scaled[j]
    #Initialization of I_0 S_0 
    S1[1, 2] = I0 #I_0
    S1[1, 1] = 1.0 - I0 #S_0
    S[1] = S1[1, 1]
    I[1] = S1[1, 2]
    h =  T^(-1)
    # SIMULATION STEP 2

    for s in 2:(T+1)
        #rand_clean = rand(Uniform(alp, 1 - alp)) # quantile(desvest[j-size_forecast:j], rand_clean)#
        desv_sim =  desvest[j] #quantile(desvest[k:j], rand_clean)
        z = rand(Normal(mean_noise, 1.0))
        #evaluate Infected next  step - rate infection with additive noise#quantile(beta[j - size_forecast:j], rand_clean)
        bt = beta[j] # quantile(beta[k:j], rand_clean)#
        S1[s, 2] = S1[s-1, 2] + (bt * h + desv_sim * z * sqrt(h)) * (1.0 - S1[s-1, 2]) * S1[s-1, 2]  - gamma * S1[s-1, 2] * h#Infected
        S1[s, 1] = 1.0 - S1[s, 2]
    end
    #saving evolution of second day (day was divided in T parts - T+1 extrems)
    I[2] = S1[(size_day + 1),2] # the forecast by day, the value 10
    S[2] = 1 - I[2]

    # Use forecasting to improve the empirical distribution
    new_I_scaled = zeros(length(I_scaled)) #historical data and forecasting
    new_I_scaled[1:j] = I_scaled[1:j]
    new_I_scaled[j+1] = I[2]#obtain new estimators of sigma and beta
    #obtain new estimators of sigma and beta
    desvest[j+1] = create_sigma_estimator(new_I_scaled[1:(j+1)],h_data, mj, (j + 1))
    beta[j+1] = create_beta_estimator(new_I_scaled[1:(j+1)],h_data, gamma, mj, (j + 1))

    #updating to recalculate next iteration
    S1 = zeros(T + 1, 2)
    S1[1, 1] = S[2]
    S1[1, 2] = I[2]
  
    for i in 3:Time_size
        for w in 2:(T+1) # quantile(beta[j - size_forecast:(j+(i-2))], rand_clean) #
            #rand_clean = rand(Uniform(alp, 1 - alp))#  quantile(desvest[j - size_forecast:(j+(i-2))], rand_clean) #
            desv_sim =  desvest[j + i - 2]#quantile(desvest[k:(j+(i-2))], rand_clean)# 
            bt =  beta[j + i - 2]# quantile(beta[k:(j+(i-2))], rand_clean)#
            z = rand(Normal(mean_noise, 1.0))
            #evaluate Infected next  step - rate infection with additive noise
            S1[w, 2] = S1[w-1, 2] + (bt* h  + desv_sim*z*sqrt(h)) * (1.0 - S1[w-1, 2]) * S1[w-1, 2] - gamma * S1[w-1, 2] * h#Infected
            #evaluate Infected next  step - rate infection with proportional noise
            S1[w, 1] = 1.0 - S1[w, 2]
        end
        #saving
        I[i] = S1[size_day + 1, 2] # save the first day, with position 10
        S[i] = 1 - I[i]

        
        new_I_scaled[j+(i-1)] = I[i]#obtain new estimators of sigma and beta
        #obtain new estimators of sigma and beta
        desvest[j+(i-1)] = create_sigma_estimator(new_I_scaled[1:(j+(i-1))],h, k, (j + (i-1) ))
        beta[j+(i-1)] = create_beta_estimator(new_I_scaled[1:(j+(i-1))],h, gamma, k, (j + (i-1)))
     
        #updating values
        S1 = zeros(T + 1, 3)
        S1[1, 1] = S[i]
        S1[1, 2] = I[i]
    end
    return [S I ]
end

#######################
# gn case
#######################

function SDE_SISg(subdivitions::Int64, size_forecast::Int64,size_day::Int64, mj::Int64, j::Int64, mean_noise::Float64, gamma::Float64, alp::Float64, I_scaled::Any)
    T = subdivitions # = 1000 
    Time_size = size_forecast # = 14
    #size_day = 5
    desvest = zeros(length(I_scaled))
    Time_interval = range(mj, j, j-mj) 
    k = mj + 1
    h_data = Time_interval[2] - Time_interval[1] 
    for i in k:j
        desvest[i] = create_sigma_estimator(I_scaled,h_data, mj, i)
    end

    beta = zeros(length(I_scaled))
    for i in k:j
        beta[i] = create_beta_estimator(I_scaled,h_data, gamma, mj, i)
    end
    S1 = zeros(T + 1, 2)#auxiliar intermediate steps: 1 (Suceptible) 2(Infected)
    #create_vectors of returns STEP 1
    S = zeros(Time_size)
    I = zeros(Time_size)
    I0 = I_scaled[j]
    #Initialization of I_0 S_0 
    S1[1, 2] = I0 #I_0
    S1[1, 1] = 1.0 - I0 #S_0
    S[1] = S1[1, 1]
    I[1] = S1[1, 2]
    h =  T^(-1)
    # SIMULATION STEP 2

    for s in 2:(T+1)
        rand_clean = rand(Uniform(alp, 1 - alp)) # desvest[j] #quantile(desvest[j-size_forecast:j], rand_clean)#
        desv_sim =  quantile(desvest[k:j], rand_clean)
        z = rand(Normal(mean_noise, 1.0))
        #evaluate Infected next  step - rate infection with additive noise# beta[j] #quantile(beta[j - size_forecast:j], rand_clean)
        bt = quantile(beta[k:j], rand_clean)#
        S1[s, 2] = S1[s-1, 2] + (bt * h + desv_sim * z * sqrt(h)) * (1.0 - S1[s-1, 2]) * S1[s-1, 2]  - gamma * S1[s-1, 2] * h#Infected
        S1[s, 1] = 1.0 - S1[s, 2]
    end
    #saving evolution of second day (day was divided in T parts - T+1 extrems)
    I[2] = S1[(size_day + 1),2] # the forecast by day, the value 10
    S[2] = 1 - I[2]

    # Use forecasting to improve the empirical distribution
    new_I_scaled = zeros(length(I_scaled)) #historical data and forecasting
    new_I_scaled[1:j] = I_scaled[1:j]
    new_I_scaled[j+1] = I[2]#obtain new estimators of sigma and beta
    #obtain new estimators of sigma and beta
    desvest[j+1] = create_sigma_estimator(new_I_scaled[1:(j+1)],h_data, mj, (j + 1))
    beta[j+1] = create_beta_estimator(new_I_scaled[1:(j+1)],h_data, gamma, mj, (j + 1))

    #updating to recalculate next iteration
    S1 = zeros(T + 1, 2)
    S1[1, 1] = S[2]
    S1[1, 2] = I[2]
  
    for i in 3:Time_size
        for w in 2:(T+1) # quantile(beta[j - size_forecast:(j+(i-2))], rand_clean) #
            rand_clean = rand(Uniform(alp, 1 - alp))#   beta[j + i - 2]#  desvest[j + i - 2]# quantile(desvest[j - size_forecast:(j+(i-2))], rand_clean) #
            desv_sim = quantile(desvest[k:(j+(i-2))], rand_clean)# 
            bt = quantile(beta[k:(j+(i-2))], rand_clean)#
            z = rand(Normal(mean_noise, 1.0))
            #evaluate Infected next  step - rate infection with additive noise
            S1[w, 2] = S1[w-1, 2] + (bt* h  + desv_sim*z*sqrt(h)) * (1.0 - S1[w-1, 2]) * S1[w-1, 2] - gamma * S1[w-1, 2] * h#Infected
            #evaluate Infected next  step - rate infection with proportional noise
            S1[w, 1] = 1.0 - S1[w, 2]
        end
        #saving
        I[i] = S1[size_day + 1, 2] # save the first day, with position 10
        S[i] = 1 - I[i]

        
        new_I_scaled[j+(i-1)] = I[i]#obtain new estimators of sigma and beta
        #obtain new estimators of sigma and beta
        desvest[j+(i-1)] = create_sigma_estimator(new_I_scaled[1:(j+(i-1))],h, k, (j + (i-1) ))
        beta[j+(i-1)] = create_beta_estimator(new_I_scaled[1:(j+(i-1))],h, gamma, k, (j + (i-1)))
     
        #updating values
        S1 = zeros(T + 1, 3)
        S1[1, 1] = S[i]
        S1[1, 2] = I[i]
    end
    return [S I ]
end

function SDE_SISg_Milstein(
    subdivitions::Int64, size_forecast::Int64, size_day::Int64,
    mj::Int64, j::Int64, mean_noise::Float64, gamma::Float64, 
    alp::Float64, I_scaled::Vector{Float64}; ϵ::Float64 = 1e-8
)
    T = subdivitions
    Time_size = size_forecast
    k = mj + 1
    h_data = (j - mj)^(-1)
    h = 1.0 / T

    desvest = zeros(length(I_scaled))
    beta = zeros(length(I_scaled))
    for i in k:j
        desvest[i] = create_sigma_estimator(I_scaled, h_data, mj, i; ϵ=ϵ)
        beta[i] = create_beta_estimator(I_scaled, h_data, gamma, mj, i)
    end

    S1 = zeros(T + 1, 2)  # 1: Susceptible, 2: Infected
    S = zeros(Time_size)
    I = zeros(Time_size)
    I0 = I_scaled[j]
    S1[1, 2] = I0
    S1[1, 1] = 1.0 - I0
    S[1] = S1[1, 1]
    I[1] = S1[1, 2]

    for s in 2:(T + 1)
        r = rand(Uniform(alp, 1 - alp))
        sigma_sim = quantile(desvest[k:j], r)
        beta_sim = quantile(beta[k:j], r)
        z = rand(Normal(mean_noise, 1.0))
        I_prev = S1[s - 1, 2]

        drift = (beta_sim * (1 - I_prev) * I_prev - gamma * I_prev)
        diffusion = sigma_sim * (1 - I_prev) * I_prev * z
        milstein_term = 0.5 * sigma_sim^2 * (1 - I_prev) * I_prev * (1 - 2 * I_prev) * (z^2 - 1)

        S1[s, 2] = I_prev + drift * h + diffusion * sqrt(h) + milstein_term * h
        S1[s, 1] = 1.0 - S1[s, 2]
    end

    I[2] = S1[size_day + 1, 2]
    S[2] = 1.0 - I[2]

    new_I_scaled = copy(I_scaled)
    new_I_scaled[j + 1] = I[2]
    desvest[j + 1] = create_sigma_estimator(new_I_scaled[1:j + 1], h_data, mj, j + 1; ϵ=ϵ)
    beta[j + 1] = create_beta_estimator(new_I_scaled[1:j + 1], h_data, gamma, mj, j + 1)

    for i in 3:Time_size
        S1[1, 1] = S[i - 1]
        S1[1, 2] = I[i - 1]
        for w in 2:(T + 1)
            r = rand(Uniform(alp, 1 - alp))
            idx_end = j + (i - 2)
            sigma_sim = quantile(desvest[k:idx_end], r)
            beta_sim = quantile(beta[k:idx_end], r)
            z = rand(Normal(mean_noise, 1.0))
            I_prev = S1[w - 1, 2]

            drift = (beta_sim * (1 - I_prev) * I_prev - gamma * I_prev)
            diffusion = sigma_sim * (1 - I_prev) * I_prev * z
            milstein_term = 0.5 * sigma_sim^2 * (1 - I_prev) * I_prev * (1 - 2 * I_prev) * (z^2 - 1)

            S1[w, 2] = I_prev + drift * h + diffusion * sqrt(h) + milstein_term * h
            S1[w, 1] = 1.0 - S1[w, 2]
        end

        I[i] = S1[size_day + 1, 2]
        S[i] = 1.0 - I[i]
        new_I_scaled[j + i - 1] = I[i]
        desvest[j + i - 1] = create_sigma_estimator(new_I_scaled[1:(j + i - 1)], h_data, k, j + i - 1; ϵ=ϵ)
        beta[j + i - 1] = create_beta_estimator(new_I_scaled[1:(j + i - 1)], h_data, gamma, k, j + i - 1)
    end

    return [S I]
end


