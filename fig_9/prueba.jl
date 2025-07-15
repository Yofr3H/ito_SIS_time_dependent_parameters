using Statistics
using Plots


# 1. Datos sintéticos con tendencia + ruido creciente
x = 1:50
y = 5 .+ 0.3 .* x .+ randn(length(x)) .* range(0.0, 2.5, length=50)  # ruido creciente

# 2. Cálculo de desviación estándar local acumulada
σ_local = [std(y[1:i]) for i in eachindex(y)]  # σ acumulativa desde el inicio
z = 1.96  # para banda de confianza del 95%

# 3. Banda creciente basada en la varianza acumulada
upper = y .+  σ_local
lower = y .-  σ_local

plotlyjs()
# 4. Graficar
Plots.plot(x, y, label="Datos observados", lw=2, legend=:topleft)
Plots.plot!(x, upper, fillrange=lower, fillalpha=0.25, label="Band CI95%", color=:green)



# 4. Graficar
Plots.plot(x, y, label="Datos observados", lw=2, legend=:topleft)
Plots.plot!(x, μ_local, label="Media local", lw=2, color=:black)
Plots.plot!(x, upper, fillrange=lower, fillalpha=0.25, label="Banda 95%", color=:blue)
