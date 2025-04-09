#calculate the forward sigma estimator for sequence of scale intected historical data I_{mj-i},...,I_{i}


function create_sigma_estimator(I_scaled::Any,h::Float64,mj::Int64, i::Int64)
    denominator = 0.0
    numerator = 0.0
    
    for k = (mj +1):i
        numerator = numerator + (I_scaled[k] - I_scaled[k - 1])^2 
        denominator = denominator + (1 - I_scaled[k-1])^2 * I_scaled[k-1]^2 * h # Δt =   1 day
    end
    return sqrt(numerator * denominator^(-1))
end
