#calculate the beta estimator for sequence of scale intected historical data I_{mj-i},...,I_{i}

function create_beta_estimator(I_scaled::Any,gamma::Float64, h::Float64, mj::Int64, i::Int64)
    # calculate beta estimator at time i
    sum_numerator = 0.0
    sum_denominator = 0.0
    
    for u in (mj):i
        sum_deltaI = 0.0
        m = 2
        for k in (mj):u
            sum_deltaI = sum_deltaI + (1 / m) * (I_scaled[k] - I_scaled[k-1])
            m = m + 1
        end
        sum_numerator = sum_numerator + sum_deltaI + gamma * I_scaled[u-1] * h # Δt =   1   day
        sum_denominator = sum_denominator + (1 - I_scaled[u-1]) * I_scaled[u-1] * h # Δt =   1   day
    end
    return sum_numerator / sum_denominator
end
