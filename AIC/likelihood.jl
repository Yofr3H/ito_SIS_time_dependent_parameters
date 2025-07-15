function conditional_log_likelihood(I::Vector{Float64}, h::Float64, gamma::Float64, mj::Int64)
    # Estimar beta y sigma primero
    T = length(I)
    beta_vals = zeros(T)
    sigma_vals = zeros(T)

    for t in (mj+1):(T - 1)
        beta_vals[t] = create_beta_estimator(I, h, gamma, mj, t)
        sigma_vals[t] = create_sigma_estimator(I, h, mj, t)
    end

    # Calcular log-verosimilitud condicional
    logL = 0.0
    for t in (mj+2):(T - 1)
        mu_t = beta_vals[t] * I[t-1] * (1 - I[t-1]) - gamma * I[t-1]
        sigma2_t = sigma_vals[t]^2
        delta_I = I[t] - I[t-1]
        logL += -0.5 * log(2 * π * sigma2_t * h) - ((delta_I - mu_t * h)^2) / (2 * sigma2_t * h)
    end

    return logL
end

function conditional_log_likelihood1(I::Vector{Float64}, h::Float64, gamma::Float64, beta::Float64, mj::Int64; ϵ::Float64 = 1e-8)
    T = length(I)
    sigma_vals = zeros(T)

    # Estimar sigma_t en cada tiempo t ≥ mj+1
    for t in (mj + 1):(T - 1)
        sigma_vals[t] = create_sigma_estimator(I, h, mj, t; ϵ=ϵ)
    end

    # Calcular log-verosimilitud condicional
    logL = 0.0
    for t in (mj + 2):(T - 1)
        I_prev = I[t - 1]
        delta_I = I[t] - I_prev
        mu_t = beta * I_prev * (1 - I_prev) - gamma * I_prev
        sigma2_t = max(sigma_vals[t]^2, ϵ)  # evitar división por cero o log(0)

        logL += -0.5 * log(2π * sigma2_t * h) - (delta_I - mu_t * h)^2 / (2 * sigma2_t * h)
    end

    return logL
end

function conditional_log_likelihood2(I::Vector{Float64}, h::Float64, gamma::Float64, mj::Int64; ϵ::Float64 = 1e-8)
    T = length(I)
    if T ≤ mj + 2
        error("Not enough data points. Need T ≥ mj + 3 to compute residuals.")
    end

    beta_vals = zeros(T)
    for t in (mj + 1):(T - 1)
        beta_vals[t] = create_beta_estimator(I, h, gamma, mj, t)
    end

    ΔI_obs = Float64[]
    ΔI_mod = Float64[]
    
    for t in (mj + 2):(T - 1)
        I_prev = I[t - 1]
        delta_obs = I[t] - I_prev
        mu_t = beta_vals[t] * I_prev * (1 - I_prev) - gamma * I_prev
        delta_mod = h * mu_t

        push!(ΔI_obs, delta_obs)
        push!(ΔI_mod, delta_mod)
    end

    residuals = ΔI_obs .- ΔI_mod
    sigma2 = var(residuals) + ϵ

    logL = -0.5 * sum(log.(2π * sigma2 * h) .+ (residuals .^ 2) ./ (sigma2 * h))

    k = T - mj
    return logL, sigma2, k
end

function conditional_log_likelihood3(I::Vector{Float64}, h::Float64, gamma::Float64, mj::Int64; ϵ::Float64 = 1e-8)
    T = length(I)
    beta_vals = zeros(T)

    # Estimar beta en cada t
    for t in (mj + 1):(T - 1)
        beta_vals[t] = create_beta_estimator(I, h, gamma, mj, t)
    end

    # Calcular mu_t y residuos
    ΔI_obs = I[(mj + 2):(T - 1)] .- I[(mj + 1):(T - 2)]
    ΔI_mod = zeros(length(ΔI_obs))
    for (i, t) in enumerate((mj + 2):(T - 1))
        I_prev = I[t - 1]
        ΔI_mod[i] = h * (beta_vals[t] * I_prev * (1 - I_prev) - gamma * I_prev)
    end

    # Estimar sigma^2 global como varianza de residuos
    residuals = ΔI_obs .- ΔI_mod
    sigma2 = var(residuals) + ϵ

    # Calcular log-verosimilitud
    logL = -0.5 * sum(log.(2π * sigma2 * h) .+ (residuals.^2) ./ (sigma2 * h))

    return logL, sigma2, T - mj  # también retorna k para usar en AIC
end





