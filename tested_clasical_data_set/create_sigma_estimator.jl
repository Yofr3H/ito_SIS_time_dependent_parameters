#calculate the forward sigma estimator for sequence of scale intected historical data I_{mj-i},...,I_{i}


function create_sigma_estimator1(I_scaled::Any,h::Float64,mj::Int64, i::Int64)
    denominator = 0.0
    numerator = 0.0
    
    for k = (mj +1):i
        numerator = numerator + (I_scaled[k] - I_scaled[k - 1])^2 
        denominator = denominator + (1 - I_scaled[k-1])^2 * I_scaled[k-1]^2 * h # Δt =   1 day
    end
    return sqrt(numerator * denominator^(-1))
end

function create_sigma_estimator(I_scaled::Vector{Float64}, h::Float64, mj::Int64, i::Int64; ϵ::Float64 = 1e-8)
    numerator = 0.0
    denominator = 0.0

    for k in (mj + 1):i
        ΔI = I_scaled[k] - I_scaled[k - 1]
        I_prev = I_scaled[k - 1]

        # Asegurar que I_prev esté en un rango válido para evitar divisiones por cero o dominio inválido
        if I_prev < ϵ || I_prev > 1 - ϵ
            continue  # Saltar este punto si es numéricamente inestable
        end

        # Estimador robusto para sigma
        numerator += ΔI^2
        denominator += ((1 - I_prev)^2 * I_prev^2 + ϵ) * h
    end

    if denominator < ϵ
        return 0.0  # retorno defensivo en caso de degeneración
    end

    return sqrt(numerator / denominator)
end
