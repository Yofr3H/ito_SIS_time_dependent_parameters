# Calculate the Root Square Mean Error RMSE between two data sequences:
# f_1, ... ,f_m  and h_1, ... ,h_m

function RSME(A::Any, B::Any)
    return sqrt(length(A))^(-1) * norm(A - B, 2)
end