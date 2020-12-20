## This whole file is outdated and not necessary after the revision of the code
# We leave it here anyway

using DrWatson
@quickactivate "PECUZAL_Julia"

using DelayEmbeddings, DynamicalSystemsBase
using DifferentialEquations
using Random
using Test
using DelimitedFiles

include("../../src/pecuzal_method.jl")

## In this script we evaluate the dependence of the returned reconstruction
# parameters on the time horizon Tw, needed for the computation of the
# L-statistic

## Dependence on p

# The time series have been computed in the script `/test/timeseries/produce_timeseries.jl`
s = readdlm("./test/timeseries/lorenz_pecora_uni_x.csv")
s = vec(s[1:5000]) # input timeseries = x component of lorenz
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
    YY, τ_valss, _, Lss , _ = pecuzal_embedding_update(s;
                                τs = 0:Tmax , w = w, samplesize = samplesize,
                                K = K, KNN = KNN, Tw = Tw)
    push!(sizes,size(YY,2))
    push!(τ_vals, τ_valss)
    push!(Ls, Lss)
end

writedlm("./scripts/computed data/dependence_on_Tw_Tws.csv",Tws)
writedlm("./scripts/computed data/dependence_on_Tw_sizes.csv",sizes)
writedlm("./scripts/computed data/dependence_on_Tw_Ls.csv",Ls)
writedlm("./scripts/computed data/dependence_on_Tw_tau_vals.csv",τ_vals)


## Plot results (see script `plot_results_params_dependencies.jl`)

# # load computed data
# Tws = readdlm("./scripts/computed data/dependence_on_Tw_Tws.csv")
# sizes = readdlm("./scripts/computed data/dependence_on_Tw_sizes.csv")
# Ls = readdlm("./scripts/computed data/dependence_on_Tw_Ls.csv")
# τ_vals = readdlm("./scripts/computed data/dependence_on_Tw_tau_vals.csv")
#
# # fill NaNs to empty entries
# τ_vals[1,3]=NaN64
# τ_vals[2,3]=NaN64
# τ_vals[3,3]=NaN64
# Ls[1,3]=NaN64
# Ls[2,3]=NaN64
# Ls[3,3]=NaN64
#
# using PyPlot
# pygui(true)
#
# fig = figure(figsize=[15,8])
#
# panelnames = ["A" "B"]
#
# lwg = 3         # linewidth of the graph
# lwa = 2         # linewidth of the axis
# fsp = 20        # Fontsize of panelnames
# fsa = 16        # Fontsize of the axis
# fsl = 12        # Fontsize of the legendentries
# fst = 16        # Fontsize of title
# axislabelsize = 16 # axislabelsize
# ticklabelsize = 12  # labelsize of ticks
# ms = 10         # markersize of maxima
# ms_style = "*"  # markerstyle
#
# subplots_adjust(left = 0.07,
# right = 0.96,
# bottom = 0.08,
# top = 0.9,
# wspace = 0.25,
# hspace = 0.25)
#
# rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
# font0 = Dict(
#         "font.size" => fsa,
#         "axes.labelweight" => "normal",
#         "axes.labelsize" => axislabelsize,
#         "axes.linewidth" => lwa,
#         "xtick.labelsize" => ticklabelsize,
#         "ytick.labelsize" => ticklabelsize,
#         "legend.fontsize" => fsl)
# merge!(rcParams, font0)
#
# ax = fig.add_subplot(1, 2, 1)
# p1 = plot(Tws, Ls[:,1], label="embedding cylce 1", linewidth=lwg)
# color1 = p1[1].get_color()
# scatter(Tws, Ls[:,1], label="")
# p2 = plot(Tws, Ls[:,2], label="embedding cylce 2", linewidth=lwg)
# color2 = p2[1].get_color()
# scatter(Tws, Ls[:,2], label="")
# legend()
# ax.text(-0.15, 1.1, panelnames[1], transform=ax.transAxes,
#      fontsize=fsp, fontweight="bold", va="top")
# title("minimum L")
# xlabel("Tw")
# ylabel(L"L_{min}")
# grid()
#
# ax = fig.add_subplot(1, 2, 2)
# plot(Tws, τ_vals[:,1], label="τ₁, embedding cylce 0", linewidth=lwg, color="r")
# scatter(Tws, τ_vals[:,1], label="", color="r")
# plot(Tws, τ_vals[:,2], label="τ₂, embedding cylce 1", linewidth=lwg, color=color1)
# scatter(Tws, τ_vals[:,2], label="", color=color1)
# plot(Tws, τ_vals[:,3], label="τ₃, embedding cylce 2", linewidth=lwg, color=color2)
# scatter(Tws, τ_vals[:,3], label="", color=color2)
# legend()
# ax.set_ylim([-1,24])
# ax.text(-0.15, 1.1, panelnames[2], transform=ax.transAxes,
#      fontsize=fsp, fontweight="bold", va="top")
# title("chosen delays")
# xlabel("Tw")
# ylabel("τ")
# grid()
