using DrWatson
@quickactivate "PECUZAL_Julia"

using DelayEmbeddings, DynamicalSystemsBase
using DifferentialEquations
using Random
using Test
using DelimitedFiles

include("../../src/pecuzal_method.jl")

## In this script we evaluate the dependence of the returned reconstruction
# parameters on the parameter KNN, needed for the computation of the
# L-statistic

## Dependence on KNN
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
K = 14
samplesize = 1
Tw = 20

Ks = 3:15
τ_vals = []
Ls = []
sizes = []
for KNN in Ks
    display(KNN)
    YY, τ_valss, _, Lss , _ = pecuzal_embedding(s;
                                τs = 0:Tmax , w = w, samplesize = samplesize,
                                K = K, KNN = KNN, Tw = Tw)
    push!(sizes,size(YY,2))
    push!(τ_vals, τ_valss)
    push!(Ls, Lss)
end

writedlm("./scripts/computed data/dependence_on_K_Ks.csv",Ks)
writedlm("./scripts/computed data/dependence_on_K_sizes.csv",sizes)
writedlm("./scripts/computed data/dependence_on_K_Ls.csv",Ls)
writedlm("./scripts/computed data/dependence_on_K_tau_vals.csv",τ_vals)


## Plot results

Ks = readdlm("./scripts/computed data/dependence_on_K_Ks.csv")
sizes = readdlm("./scripts/computed data/dependence_on_K_sizes.csv")
Ls = readdlm("./scripts/computed data/dependence_on_K_Ls.csv")
τ_vals = readdlm("./scripts/computed data/dependence_on_K_tau_vals.csv")

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
p1 = plot(Ks, Ls[:,1], label="embedding cylce 1", linewidth=lwg)
color1 = p1[1].get_color()
scatter(Ks, Ls[:,1], label="")
p2 = plot(Ks, Ls[:,2], label="embedding cylce 2", linewidth=lwg)
color2 = p2[1].get_color()
scatter(Ks, Ls[:,2], label="")
legend()
ax.text(-0.15, 1.1, panelnames[1], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("minimum L")
xlabel("K-NN")
ylabel(L"L_{min}")
grid()

ax = fig.add_subplot(1, 2, 2)
plot(Ks, τ_vals[:,1], label="τ₁, embedding cylce 0", linewidth=lwg, color="r")
scatter(Ks, τ_vals[:,1], label="", color="r")
plot(Ks, τ_vals[:,2], label="τ₂, embedding cylce 1", linewidth=lwg, color=color1)
scatter(Ks, τ_vals[:,2], label="", color=color1)
plot(Ks, τ_vals[:,3], label="τ₃, embedding cylce 2", linewidth=lwg, color=color2)
scatter(Ks, τ_vals[:,3], label="", color=color2)
legend()
ax.set_ylim([-1,24])
ax.text(-0.15, 1.1, panelnames[2], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("chosen delays")
xlabel("k-NN")
ylabel("τ")
grid()
