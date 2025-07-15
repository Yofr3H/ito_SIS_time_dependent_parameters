include("./likelihood.jl")

function compute_AIC(I::Vector{Float64}, h::Float64, gamma::Float64, mj::Int64)
    logL = conditional_log_likelihood(I, h, gamma, mj)

    # Number of parameters:  beta and  sigma for each time t ∈ [mj+1, T-1]
    T = length(I)
    k =  2 * (T - mj - 1)  # beta_t and sigma_t by time #k=0 case 

    AIC = 2k - 2 * logL
    return AIC
end

function compute_AIC1(I::Vector{Float64}, h::Float64, gamma::Float64, beta::Float64, mj::Int64)
    logL = conditional_log_likelihood1(I, h, gamma,beta,mj)

    # Number of parameters:  beta and  sigma for each time t ∈ [mj+1, T-1]
    T = length(I)
    k = (T - mj - 2)  # number of estimated sigma_t  

    AIC = 2k - 2 * logL
    return AIC
end


function compute_AIC2(I::Vector{Float64}, h::Float64, gamma::Float64, mj::Int64)
    logL, sigma2, k = conditional_log_likelihood3(I, h, gamma, mj)
    AIC = 2 * k - 2 * logL
    return AIC
end


function AIC_SIS_deterministic(
    I_data::Vector{Float64}, mj::Int, steps::Int;
    beta::Float64 = 0.10, gamma::Float64 = 1/55,
    h::Float64 = 1.0, sigma2::Union{Nothing,Float64} = nothing
    )
    # simulation trajectory
    I_model = SIS_deterministic_forecast(I_data, mj, steps; beta=beta, gamma=gamma, h=h)

    # data
    I_obs = I_data[mj:(mj + steps - 1)]

    # increments
    ΔI_obs = diff(I_obs)
    ΔI_mod = diff(I_model)

    #  variance estimation
    if sigma2 === nothing
        residuals = ΔI_obs .- ΔI_mod
        sigma2 = var(residuals) + 1e-8  # prevent degeneration
    end

    # Log-verosimilitud under condicional gaussian assumption
    logL = -0.5 * sum(log.(2π * sigma2 * h) .+ ((ΔI_obs .- ΔI_mod).^2) ./ (sigma2 * h))

    # Number of parameters estimated (only beta, gamma)
    k = 2

    # AIC
    AIC = 2k - 2 * logL
    return AIC
end


