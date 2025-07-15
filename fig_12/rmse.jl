# Calculate the Root Square Mean Error RMSE between two data sequences:
# f_1, ... ,f_m  and h_1, ... ,h_m

function RSME(A::Any, B::Any)
    return sqrt(length(A))^(-1) * norm(A - B, 2)
end

function MAPE(y_true::Vector{Float64}, y_pred::Vector{Float64})
    if length(y_true) != length(y_pred)
        error("Vectors must have the same length")
    end
    return 100 * mean(abs.((y_true .- y_pred) ./ y_true))
end