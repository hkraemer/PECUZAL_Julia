using DrWatson
@quickactivate "PECUZAL_Julia"

using DelayEmbeddings, DynamicalSystemsBase
using DifferentialEquations
using Random
using Test
using DelimitedFiles

include("../../src/pecuzal_method.jl")

## In this script we evaluate the dependence of the returned reconstruction
# parameters on the parameter α, needed for the computation of the
# contiuity-statistic

## Dependence on α

# The time series have been computed in the script `/test/timeseries/produce_timeseries.jl`
s = readdlm("./test/timeseries/lorenz_pecora_uni_x.csv")
s = vec(s[1:5000]) # input timeseries = x component of lorenz
w = estimate_delay(s, "mi_min")
Tmax = 100
KNN = 3
samplesize = 1
Tw = 20
K = 14
alphas = [0.05 0.01]

τ_vals = []
Ls = []
sizes = []
εs = []
for alpha in alphas
    display(alpha)
    YY, τ_valss, _, Lss , ε = pecuzal_embedding_update(s;
                                τs = 0:Tmax , w = w, samplesize = samplesize,
                                K = K, KNN = KNN, Tw = Tw, α = alpha)
    push!(sizes,size(YY,2))
    push!(τ_vals, τ_valss)
    push!(Ls, Lss)
    push!(εs, ε)
end

writedlm("./scripts/computed data/dependence_on_alpha_alphas.csv",alphas)
writedlm("./scripts/computed data/dependence_on_alpha_sizes.csv",sizes)
writedlm("./scripts/computed data/dependence_on_alpha_Ls.csv",Ls)
writedlm("./scripts/computed data/dependence_on_alpha_tau_vals.csv",τ_vals)
writedlm("./scripts/computed data/dependence_on_alpha_epsilons.csv",εs)

## Plot results (see script `plot_results_params_dependencies.jl`)

# alphas = readdlm("./scripts/computed data/dependence_on_alpha_alphas.csv")
# sizes = readdlm("./scripts/computed data/dependence_on_alpha_sizes.csv")
# Ls = readdlm("./scripts/computed data/dependence_on_alpha_Ls.csv")
# τ_vals = readdlm("./scripts/computed data/dependence_on_alpha_tau_vals.csv")
# εs = readdlm("./scripts/computed data/dependence_on_alpha_epsilons.csv")
#
# using PyPlot
# pygui(true)
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
# taus = vec(0:Tmax)
# fig = figure(figsize=[15,8])
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
# ylims = [0, 1.4]
#
# ax = fig.add_subplot(1, 2, 1)
# plot(taus,εs[1][:,1], label="τs=$(τ_vals[1][1])", linewidth = lwg)
# plot(taus,εs[1][:,2], label="τs=$(τ_vals[1][2])", linewidth = lwg)
# plot(taus,εs[1][:,3], label="τs=$(τ_vals[1][3])", linewidth = lwg)
# plot([taus[τ_vals[1][2]+1]], [εs[1][τ_vals[1][2]+1,1]], ms_style, label="", markersize = ms, color = "black")
# plot([taus[τ_vals[1][3]+1]], [εs[1][τ_vals[1][3]+1,2]], ms_style, label="", markersize = ms, color = "black")
# legend()
# ax.set_ylim(ylims)
# ax.text(-0.15, 1.1, panelnames[1], transform=ax.transAxes,
#      fontsize=fsp, fontweight="bold", va="top")
# title("α=0.05")
# xlabel("τ")
# ylabel("⟨ε★⟩")
# grid()
#
# ax = fig.add_subplot(1, 2, 2)
# plot(taus,εs[2][:,1], label="τs=$(τ_vals[2][1])", linewidth = lwg)
# plot(taus,εs[2][:,2], label="τs=$(τ_vals[2][2])", linewidth = lwg)
# plot(taus,εs[2][:,3], label="τs=$(τ_vals[2][3])", linewidth = lwg)
# plot([taus[τ_vals[2][2]+1]], [εs[2][τ_vals[2][2]+1,1]], ms_style, label="", markersize = ms, color = "black")
# plot([taus[τ_vals[2][3]+1]], [εs[2][τ_vals[2][3]+1,2]], ms_style, label="", markersize = ms, color = "black")
# legend()
# ax.set_ylim(ylims)
# ax.text(-0.15, 1.1, panelnames[2], transform=ax.transAxes,
#      fontsize=fsp, fontweight="bold", va="top")
# title("α=0.01")
# xlabel("τ")
# grid()
