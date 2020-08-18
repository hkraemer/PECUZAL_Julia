using DrWatson
@quickactivate "new-embedding-methods"

using DelayEmbeddings, DynamicalSystemsBase
using DifferentialEquations
using Random
using Test

include("/src/pecora_uzal_method.jl")

Random.seed!(414515)
lor = Systems.lorenz([2*rand(1); 2*rand(1); 2*rand(1)]; ρ=60)
data = trajectory(lor, 50; dt=0.01, Ttr = 10)

s = data[:, 1] # input timeseries = y component of lorenz
w = estimate_delay(s, "mi_min")
Tmax = 50
K = 14
samplesize = 1
KNN = 3


@time Y, τ_vals, ts_vals, Ls , ε★ = pecora_uzal_embedding(s;
                                    τs = 0:Tmax , w = w, samplesize = samplesize,
                                    K = K, KNN = KNN)

using PyPlot
pygui(true)
figure()
for i in 1:length(ts_vals)
    plot(ε★[:, i], label = "τs = $(τ_vals[1:i])")
end
ylabel("⟨ε★⟩")
title("lorenz system")
legend()
grid()
