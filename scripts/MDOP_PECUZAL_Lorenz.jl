using DrWatson
@quickactivate "new-embedding-methods"

using DelayEmbeddings, DynamicalSystemsBase
using DifferentialEquations
using Random
using Test

include("/Users/hkraemer/Documents/Git/new-embedding-methods/new-embedding-methods/src/pecora_uzal_method.jl")

## In this script we compare the reconstruction results gained from MDOP- and from
# Pecora-Uzal-embedding for the univariate Lorenz system

## Test MDOP and Pecora-Uzal on univariate Lorenz

Random.seed!(414515)
lor = Systems.lorenz([2*rand(1); 2*rand(1); 2*rand(1)]; ρ=60)
data = trajectory(lor, 50; dt=0.01, Ttr = 10)

s = data[:, 1] # input timeseries = x component of lorenz
w = estimate_delay(s, "mi_min")
Tmax = 100
K = 14
samplesize = 1
KNN = 3
Tw = 50


@time Y, τ_vals, ts_vals, Ls , ε★ = pecora_uzal_embedding(s;
                                    τs = 0:Tmax , w = w, samplesize = samplesize,
                                    K = K, KNN = KNN, Tw = Tw)

# TODO make mdop_maximum_delay work
# @time Tmax2 = mdop_maximum_delay(s)

Tmax2 = 20
@time Y2, τ_vals2, ts_vals2, FNNs , βs = mdop_embedding(s;
                                    τs = 0:Tmax2 , w = w)

# Plot results
using PyPlot
pygui(true)
figure(figsize=[12, 8])
subplot(121)
for i in 1:length(ts_vals)
    plot(ε★[:, i], label = "τs = $(τ_vals[1:i])")
end
ylabel("⟨ε★⟩")
title("Pecora-Uzal (Lorenz system)")
legend()
grid()

subplot(122)
for i in 1:length(ts_vals2)
    plot(βs[:, i], label = "τs = $(τ_vals2[1:i])")
end
ylabel("β-statistic")
title("MDOP (Lorenz system)")
legend()
grid()

fig = figure(figsize=[20, 10])
ax = fig.add_subplot(1, 2, 1, projection="3d")
ax.plot3D(Y[:,1], Y[:,2], Y[:,3])
xlabel("x(t+$(τ_vals[1]))")
ylabel("x(t+$(τ_vals[2]))")
zlabel("x(t+$(τ_vals[3]))")
title("Pecora-Uzal Attractor (Lorenz system) \n τs=$(τ_vals), ts=$(ts_vals)")
grid()

ax = fig.add_subplot(1, 2, 2, projection="3d")
ax.plot3D(Y2[:,1], Y2[:,2], Y2[:,3])
xlabel("x(t+$(τ_vals2[1]))")
ylabel("x(t+$(τ_vals2[2]))")
zlabel("x(t+$(τ_vals2[3]))")
title("MDOP Attractor (Lorenz system) \n τs=$(τ_vals2), ts=$(ts_vals2)")
grid()
