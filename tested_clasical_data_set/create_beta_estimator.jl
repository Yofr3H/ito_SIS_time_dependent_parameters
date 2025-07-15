#calculate the beta estimator for sequence of scale intected historical data I_{mj-i},...,I_{i}

function create_beta_estimator1(I::Any,h::Float64,gamma::Float64, mj::Int64, i::Int64)
    # calculate sigma estimator at time i
    beta_estimator = 0.0
      
    for j in (mj + 1):i
        beta_estimator = beta_estimator +  (I[j] - I[j-1]) * ( I[j - 1] * (1 - I[j-1] ) )^(-1) + gamma * h * ( (1 - I[j - 1] ) )^(-1)
    end
    return beta_estimator * (h*(i-mj))
end

function create_beta_estimator(I::AbstractVector{Float64}, h::Float64, gamma::Float64, mj::Int64, i::Int64)
    beta_estimator = 0.0
    eps = 1e-8  # Umbral para evitar divisiones por cero

    for j in (mj + 1):i
        I_prev = I[j - 1]
        denominator = max(I_prev * (1.0 - I_prev), eps)
        correction = gamma * h / max(1.0 - I_prev, eps)

        beta_estimator += (I[j] - I_prev) / denominator + correction
    end

    interval_length = max(h * (i - mj), eps)  # Evita división por cero si el intervalo es nulo
    return beta_estimator / interval_length
end
