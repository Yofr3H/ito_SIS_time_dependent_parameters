using Statistics
using Distributions

function IC(y::Any)
 
    σ_local = [std(y[1:i]) for i in eachindex(y)]
    σ_local[1] = 0.0 
    upper = y .+  σ_local
    lower = y .-  σ_local
    return [ upper, lower]
end
