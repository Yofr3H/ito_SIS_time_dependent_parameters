#calculate the forward sigma estimator for sequence of scale intected historical data I_{mj-i},...,I_{i}


function create_sigma_estimator1(I_scaled::Any, h::Float64, mj::Int64, i::Int64)
    denominator = 0.0
    numerator = 0.0
    
    for k = mj:i
        sum_deltaI = 0.0
        m = 2
        for s = mj:k
            sum_deltaI = sum_deltaI +  (I_scaled[s] - I_scaled[s-1])^2 # Δt =   1 day (1 / m) *
            m = m + 1
        end
        numerator = numerator + sum_deltaI
        denominator = denominator + (1 - I_scaled[k])^2 * I_scaled[k]^2 * h # Δt =   1 day
    end
    return sqrt(numerator / denominator)
end

function create_sigma_estimator(I_scaled::Any, i::Int64)
    denominator = 0.0
    numerator = 0.0
    
    for k = 2:i
        numerator = numerator +  (I_scaled[k] - I_scaled[k-1])^2 # Δt =   1 day (1 / m) *
      
        denominator = denominator + (1 - I_scaled[k])^2 * I_scaled[k]^2 + (1 - I_scaled[k-1])^2 * I_scaled[k-1]^2 # Δt =   1 day
    end
    return sqrt(numerator / denominator)
    #return sqrt(numerator)
end
