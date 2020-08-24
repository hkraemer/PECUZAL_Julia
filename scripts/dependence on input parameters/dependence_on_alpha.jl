using DrWatson
@quickactivate "new-embedding-methods"

using DelayEmbeddings, DynamicalSystemsBase
using DifferentialEquations
using Random
using Test
using DelimitedFiles

include("/Users/hkraemer/Documents/Git/new-embedding-methods/new-embedding-methods/src/pecora_uzal_method.jl")

## In this script we evaluate the dependence of the returned reconstruction
# parameters on the time horizon Tw, needed for the computation of the
# L-statistic

## Dependence on Tw
# For comparison reasons using Travis CI we carry out the integration on a UNIX
# OS and save the resulting time series
# lo = Systems.lorenz([1.0, 1.0, 50.0])
# tr = trajectory(lo, 100; dt = 0.01, Ttr = 10)
# x = tr[:, 1]
# writedlm("lorenz_uzal_2.csv", x)

s = readdlm("lorenz_uzal_2.csv")
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
    YY, τ_valss, _, Lss , ε = pecora_uzal_embedding(s;
                                τs = 0:Tmax , w = w, samplesize = samplesize,
                                K = K, KNN = KNN, Tw = Tw, α = alpha)
    push!(sizes,size(YY,2))
    push!(τ_vals, τ_valss)
    push!(Ls, Lss)
    push!(εs, ε)
end

writedlm("dependence_on_alpha_alphas.csv",alphas)
writedlm("dependence_on_alpha_sizes.csv",sizes)
writedlm("dependence_on_alpha_Ls.csv",Ls)
writedlm("dependence_on_alpha_tau_vals.csv",τ_vals)
writedlm("dependence_on_alpha_epsilons.csv",εs)

## Plot results

alphas = readdlm("/Users/hkraemer/Documents/Git/new-embedding-methods/new-embedding-methods/scripts/computed data/dependence_on_alpha_alphas.csv")
sizes = readdlm("/Users/hkraemer/Documents/Git/new-embedding-methods/new-embedding-methods/scripts/computed data/dependence_on_alpha_sizes.csv")
Ls = readdlm("/Users/hkraemer/Documents/Git/new-embedding-methods/new-embedding-methods/scripts/computed data/dependence_on_alpha_Ls.csv")
τ_vals = readdlm("/Users/hkraemer/Documents/Git/new-embedding-methods/new-embedding-methods/scripts/computed data/dependence_on_alpha_tau_vals.csv")
εs = readdlm("/Users/hkraemer/Documents/Git/new-embedding-methods/new-embedding-methods/scripts/computed data/dependence_on_alpha_epsilons.csv")

using PyPlot
pygui(true)

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

taus = vec(0:Tmax)
fig = figure(figsize=[15,8])

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

ylims = [0, 1.4]

ax = fig.add_subplot(1, 2, 1)
plot(taus,εs[1][:,1], label="τs=$(τ_vals[1][1])", linewidth = lwg)
plot(taus,εs[1][:,2], label="τs=$(τ_vals[1][2])", linewidth = lwg)
plot(taus,εs[1][:,3], label="τs=$(τ_vals[1][3])", linewidth = lwg)
plot([taus[τ_vals[1][2]+1]], [εs[1][τ_vals[1][2]+1,1]], ms_style, label="", markersize = ms, color = "black")
plot([taus[τ_vals[1][3]+1]], [εs[1][τ_vals[1][3]+1,2]], ms_style, label="", markersize = ms, color = "black")
legend()
ax.set_ylim(ylims)
ax.text(-0.15, 1.1, panelnames[1], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("α=0.05")
xlabel("τ")
ylabel("⟨ε★⟩")
grid()

ax = fig.add_subplot(1, 2, 2)
plot(taus,εs[2][:,1], label="τs=$(τ_vals[2][1])", linewidth = lwg)
plot(taus,εs[2][:,2], label="τs=$(τ_vals[2][2])", linewidth = lwg)
plot(taus,εs[2][:,3], label="τs=$(τ_vals[2][3])", linewidth = lwg)
plot([taus[τ_vals[2][2]+1]], [εs[2][τ_vals[2][2]+1,1]], ms_style, label="", markersize = ms, color = "black")
plot([taus[τ_vals[2][3]+1]], [εs[2][τ_vals[2][3]+1,2]], ms_style, label="", markersize = ms, color = "black")
legend()
ax.set_ylim(ylims)
ax.text(-0.15, 1.1, panelnames[2], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("α=0.01")
xlabel("τ")
grid()
