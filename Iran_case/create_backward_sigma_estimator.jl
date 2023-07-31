
function create_backward_sigma_estimator(I_scaled::Any, mj::Int64, i::Int64)
    denominator = 0.0
    numerator = 0.0
    for k = i:(-1):(mj) # (i-mj+1):i
        sum_deltaI = 0.0
        m = 2
        for s = k:(-1):(mj) #
            sum_deltaI = sum_deltaI + (1 / m) * (I_scaled[s] - I_scaled[s-1])^2 # Δt =   1 day
            m = m + 1
        end
        numerator = numerator + sum_deltaI
        denominator = denominator + (1 - I_scaled[k])^2 * I_scaled[k]^2 # Δt =   1 day
    end
    return sqrt(numerator / denominator)
end
