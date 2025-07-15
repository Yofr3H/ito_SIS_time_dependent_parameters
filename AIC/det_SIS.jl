function SIS_deterministic_forecast(I_data::Vector{Float64}, mj::Int, steps::Int; 
    beta::Float64 = 0.10, gamma::Float64 = 1/55, h::Float64 = 1.0)
    # to save trajectory
    I = zeros(steps)
    
    # initial value - historical data
    I[1] = I_data[mj]

    #  ODE SIS deterministic -  Euler
    for t in 2:steps
        It = I[t - 1]
        dI = (beta * It * (1 - It) - gamma * It) * h
        I[t] = It + dI
        I[t] = clamp(I[t], 0.0, 1.0)  # limitar a [0,1]
    end
    
    return I
end

