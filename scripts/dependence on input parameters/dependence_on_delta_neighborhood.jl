using DrWatson
@quickactivate "new-embedding-methods"

using DelayEmbeddings, DynamicalSystemsBase
using DifferentialEquations
using Random
using Test
using DelimitedFiles

include("../../src/pecora_uzal_method.jl")

## In this script we evaluate the dependence of the returned reconstruction
# parameters on the parameter determining the δ-neighborhood, needed for the
# computation of the continuity-statistic

## Dependence on δ-neighborhood
# For comparison reasons using Travis CI we carry out the integration on a UNIX
# OS and save the resulting time series
# lo = Systems.lorenz([1.0, 1.0, 50.0])
# tr = trajectory(lo, 100; dt = 0.01, Ttr = 10)
# x = tr[:, 1] # x-component of time series
# y = tr[:, 2] #y-component of time series
# writedlm("./test/timeseries/lorenz_pecora_uni_x.csv", x)
# writedlm("./test/timeseries/lorenz_pecora_uni_y.csv", y)
# writedlm("./test/timeseries/lorenz_pecora_multi.csv", tr)

s = readdlm("./test/timeseries/lorenz_pecora_uni_x.csv")
s = vec(s[1:5000]) # input timeseries = x component of lorenz
w = estimate_delay(s, "mi_min")
Tmax = 100
KNN = 3
samplesize = 1
Tw = 20

Ks = 8:20
τ_vals = []
Ls = []
sizes = []
for K in Ks
    display(K)
    YY, τ_valss, _, Lss , _ = pecuzal_embedding(s;
                                τs = 0:Tmax , w = w, samplesize = samplesize,
                                K = K, KNN = KNN, Tw = Tw)
    push!(sizes,size(YY,2))
    push!(τ_vals, τ_valss)
    push!(Ls, Lss)
end

writedlm("./scripts/computed data/dependence_on_delta_neighborhood_Ks.csv",Ks)
writedlm("./scripts/computed data/dependence_on_delta_neighborhood_sizes.csv",sizes)
writedlm("./scripts/computed data/dependence_on_delta_neighborhood_Ls.csv",Ls)
writedlm("./scripts/computed data/dependence_on_delta_neighborhood_tau_vals.csv",τ_vals)

## Plot results

deltas = readdlm("./scripts/computed data/dependence_on_delta_neighborhood_Ks.csv")
sizes = readdlm("./scripts/computed data/dependence_on_delta_neighborhood_sizes.csv")
Ls = readdlm("./scripts/computed data/dependence_on_delta_neighborhood_Ls.csv")
τ_vals = readdlm("./scripts/computed data/dependence_on_delta_neighborhood_tau_vals.csv")


using PyPlot
pygui(true)

fig = figure(figsize=[15,8])

panelnames = ["A" "B"]

lwg = 3         # linewidth of the graph
lwa = 2         # linewidth of the axis
fsp = 20        # Fontsize of panelnames
fsa = 16        # Fontsize of the axis
fsl = 12        # Fontsize of the legendentries
fst = 16        # Fontsize of title
axislabelsize = 16 # axislabelsize
ticklabelsize = 12  # labelsize of ticks
ms = 10         # markersize of maxima
ms_style = "*"  # markerstyle

subplots_adjust(left = 0.07,
right = 0.96,
bottom = 0.08,
top = 0.9,
wspace = 0.25,
hspace = 0.25)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
font0 = Dict(
        "font.size" => fsa,
        "axes.labelweight" => "normal",
        "axes.labelsize" => axislabelsize,
        "axes.linewidth" => lwa,
        "xtick.labelsize" => ticklabelsize,
        "ytick.labelsize" => ticklabelsize,
        "legend.fontsize" => fsl)
merge!(rcParams, font0)

ax = fig.add_subplot(1, 2, 1)
p1 = plot(deltas, Ls[:,1], label="embedding cylce 1", linewidth=lwg)
color1 = p1[1].get_color()
scatter(deltas, Ls[:,1], label="")
p2 = plot(deltas, Ls[:,2], label="embedding cylce 2", linewidth=lwg)
color2 = p2[1].get_color()
scatter(deltas, Ls[:,2], label="")
legend()
ax.text(-0.15, 1.1, panelnames[1], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("minimum L")
xlabel(L"\delta-Neighborhood size")
ylabel(L"L_{min}")
grid()

ax = fig.add_subplot(1, 2, 2)
plot(deltas, τ_vals[:,1], label="τ₁, embedding cylce 0", linewidth=lwg, color="r")
scatter(deltas, τ_vals[:,1], label="", color="r")
plot(deltas, τ_vals[:,2], label="τ₂, embedding cylce 1", linewidth=lwg, color=color1)
scatter(deltas, τ_vals[:,2], label="", color=color1)
plot(deltas, τ_vals[:,3], label="τ₃, embedding cylce 2", linewidth=lwg, color=color2)
scatter(deltas, τ_vals[:,3], label="", color=color2)
legend()
ax.set_ylim([-1,24])
ax.text(-0.15, 1.1, panelnames[2], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("chosen delays")
xlabel(L"\delta-Neighborhood size")
ylabel("τ")
grid()
