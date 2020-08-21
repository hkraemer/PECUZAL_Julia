using DrWatson
@quickactivate "new-embedding-methods"

using DelayEmbeddings, DynamicalSystemsBase
using DifferentialEquations
using Random
using Test

include("/Users/hkraemer/Documents/Git/new-embedding-methods/new-embedding-methods/src/pecora_uzal_method.jl")

## In this script we evaluate the dependence of the returned reconstruction
# parameters on the time horizon Tw, needed for the computation of the
# L-statistic

## Dependence on Tw
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

Tws = 1:100
τ_vals = []
Ls = []
sizes = []
for Tw in Tws
    display(Tw)
    YY, τ_valss, _, Lss , _ = pecora_uzal_embedding(s;
                                τs = 0:Tmax , w = w, samplesize = samplesize,
                                K = K, KNN = KNN, Tw = Tw)
    push!(sizes,size(YY,2))
    push!(τ_vals, τ_valss)
    push!(Ls, Lss)
end

# Display results
using PyPlot
pygui(true)

fig = figure(figsize=[15,8])
ax = fig.add_subplot(1, 2, 1)
plot(Tws,sizes)
ax.set_xticks(Tws[1:10:end])
ax.set_xticklabels(Tws[1:10:end])
title("encountered embedding dimension")
xlabel("Tw")
ylabel("m")
grid()

ax = fig.add_subplot(1, 2, 2)
ax.plot(Tws, minimum.(Ls))
ax.set_xticks(Tws[1:10:end])
ax.set_xticklabels(Tws[1:10:end])
title("minimum L-vals")
xlabel("Tw")
ylabel("L")
grid()

width = 0.3

fig = figure(figsize=[20,10])
ax = fig.add_subplot(111)
for i = 2:length(Tws)
ax.bar(Tws[i],τ_vals[i][1], width, color="r")
ax.bar(Tws[i]+width,τ_vals[i][2], width, color="g")
ax.bar(Tws[i]+width*2,τ_vals[i][3], width, color="b")
end
# ax.set_xticks(Tws[1:10:end])
# ax.set_xticklabels(Tws[1:10:end])
ax.set_ylim((-1,14))
title("τ-vals")
xlabel("Tw")
ylabel("τ")
grid()
