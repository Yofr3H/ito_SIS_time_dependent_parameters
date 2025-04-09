#calculate the beta estimator for sequence of scale intected historical data I_{mj-i},...,I_{i}

function create_beta_estimator(I::Any,h::Float64,gamma::Float64, mj::Int64, i::Int64)
    # calculate sigma estimator at time i
    beta_estimator = 0.0
      
    for j in (mj + 1):i
        beta_estimator = beta_estimator +  (I[j] - I[j-1]) * ( I[j - 1] * (1 - I[j-1] ) )^(-1) + gamma * h * ( (1 - I[j - 1] ) )^(-1)
    end
    return beta_estimator * (h*(i-mj))
end
