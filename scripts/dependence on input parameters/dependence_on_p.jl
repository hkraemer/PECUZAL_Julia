using DrWatson
@quickactivate "new-embedding-methods"

using DelayEmbeddings, DynamicalSystemsBase
using DifferentialEquations
using Random
using Test
using DelimitedFiles

include("../../src/pecora_uzal_method.jl")

## In this script we evaluate the dependence of the returned reconstruction
# parameters on the binomial parameter p, needed for the computation of the
# continuity-statistic

## Dependence on p
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
K = 14
ps = [0.6 0.5 0.4 0.3 0.2 0.1]

τ_vals = []
Ls = []
sizes = []
εs = []
for p in ps
    display(p)
    YY, τ_valss, _, Lss , ε = pecuzal_embedding(s;
                                τs = 0:Tmax , w = w, samplesize = samplesize,
                                K = K, KNN = KNN, Tw = Tw, p = p)
    push!(sizes,size(YY,2))
    push!(τ_vals, τ_valss)
    push!(Ls, Lss)
    push!(εs, ε)
end

writedlm("./scripts/computed data/dependence_on_p_ps.csv",ps)
writedlm("./scripts/computed data/dependence_on_p_sizes.csv",sizes)
writedlm("./scripts/computed data/dependence_on_p_Ls.csv",Ls)
writedlm("./scripts/computed data/dependence_on_p_tau_vals.csv",τ_vals)
writedlm("./scripts/computed data/dependence_on_p_epsilons.csv",εs)

## Plot results

ps = readdlm("./scripts/computed data/dependence_on_p_ps.csv")
sizes = readdlm("./scripts/computed data/dependence_on_p_sizes.csv")
Ls = readdlm("./scripts/computed data/dependence_on_p_Ls.csv"
τ_vals = readdlm("./scripts/computed data/dependence_on_p_tau_vals.csv")
εs = readdlm("./scripts/computed data/dependence_on_p_epsilons.csv")

using PyPlot
pygui(true)

taus = vec(0:Tmax)

panelnames = ["A" "B" "C" "D" "E" "F"]

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

fig = figure(figsize=[15,8])
subplots_adjust(left = 0.08,
right = 0.97,
bottom = 0.07,
top = 0.93,
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

ax = fig.add_subplot(2, 3, 1)
plot(taus,εs[1][:,1], label="τs=$(τ_vals[1][1])", linewidth=lwg)
plot(taus,εs[1][:,2], label="τs=$(τ_vals[1][2])", linewidth=lwg)
plot(taus,εs[1][:,3], label="τs=$(τ_vals[1][3])", linewidth=lwg)
plot([taus[τ_vals[1][2]+1]], [εs[1][τ_vals[1][2]+1,1]], ms_style, label="", markersize = ms, color = "black")
plot([taus[τ_vals[1][3]+1]], [εs[1][τ_vals[1][3]+1,2]], ms_style, label="", markersize = ms, color = "black")
legend()
ax.set_ylim(ylims)
ax.text(-0.15, 1.15, panelnames[1], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("p=0.6")
ylabel("⟨ε★⟩")
grid()

ax = fig.add_subplot(2, 3, 2)
plot(taus,εs[2][:,1], label="τs=$(τ_vals[2][1])", linewidth=lwg)
plot(taus,εs[2][:,2], label="τs=$(τ_vals[2][2])", linewidth=lwg)
plot(taus,εs[2][:,3], label="τs=$(τ_vals[2][3])", linewidth=lwg)
plot([taus[τ_vals[2][2]+1]], [εs[2][τ_vals[2][2]+1,1]], ms_style, label="", markersize = ms, color = "black")
plot([taus[τ_vals[2][3]+1]], [εs[2][τ_vals[2][3]+1,2]], ms_style, label="", markersize = ms, color = "black")
legend()
ax.set_ylim(ylims)
ax.text(-0.15, 1.15, panelnames[2], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("p=0.5")
grid()

ax = fig.add_subplot(2, 3, 3)
plot(taus,εs[3][:,1], label="τs=$(τ_vals[3][1])", linewidth=lwg)
plot(taus,εs[3][:,2], label="τs=$(τ_vals[3][2])", linewidth=lwg)
plot(taus,εs[3][:,3], label="τs=$(τ_vals[3][3])", linewidth=lwg)
plot([taus[τ_vals[3][2]+1]], [εs[3][τ_vals[3][2]+1,1]], ms_style, label="", markersize = ms, color = "black")
plot([taus[τ_vals[3][3]+1]], [εs[3][τ_vals[3][3]+1,2]], ms_style, label="", markersize = ms, color = "black")
legend()
ax.set_ylim(ylims)
ax.text(-0.15, 1.15, panelnames[3], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("p=0.4")
grid()

ax = fig.add_subplot(2, 3, 4)
plot(taus,εs[4][:,1], label="τs=$(τ_vals[4][1])", linewidth=lwg)
plot(taus,εs[4][:,2], label="τs=$(τ_vals[4][2])", linewidth=lwg)
plot(taus,εs[4][:,3], label="τs=$(τ_vals[4][3])", linewidth=lwg)
plot([taus[τ_vals[4][2]+1]], [εs[4][τ_vals[4][2]+1,1]], ms_style, label="", markersize = ms, color = "black")
plot([taus[τ_vals[4][3]+1]], [εs[4][τ_vals[4][3]+1,2]], ms_style, label="", markersize = ms, color = "black")
legend()
ax.set_ylim(ylims)
ax.text(-0.15, 1.15, panelnames[4], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("p=0.3")
xlabel("τ")
ylabel("⟨ε★⟩")
grid()

ax = fig.add_subplot(2, 3, 5)
plot(taus,εs[5][:,1], label="τs=$(τ_vals[5][1])", linewidth=lwg)
plot(taus,εs[5][:,2], label="τs=$(τ_vals[5][2])", linewidth=lwg)
plot(taus,εs[5][:,3], label="τs=$(τ_vals[5][3])", linewidth=lwg)
plot([taus[τ_vals[5][2]+1]], [εs[5][τ_vals[5][2]+1,1]], ms_style, label="", markersize = ms, color = "black")
plot([taus[τ_vals[5][3]+1]], [εs[5][τ_vals[5][3]+1,2]], ms_style, label="", markersize = ms, color = "black")
legend()
ax.set_ylim(ylims)
ax.text(-0.15, 1.15, panelnames[5], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("p=0.2")
xlabel("τ")
grid()

ax = fig.add_subplot(2, 3, 6)
plot(taus,εs[6][:,1], label="τs=$(τ_vals[6][1])", linewidth=lwg)
plot(taus,εs[6][:,2], label="τs=$(τ_vals[6][2])", linewidth=lwg)
plot(taus,εs[6][:,3], label="τs=$(τ_vals[6][3])", linewidth=lwg)
plot([taus[τ_vals[6][2]+1]], [εs[6][τ_vals[6][2]+1,1]], ms_style, label="", markersize = ms, color = "black")
plot([taus[τ_vals[6][3]+1]], [εs[6][τ_vals[6][3]+1,2]], ms_style, label="", markersize = ms, color = "black")
legend()
ax.set_ylim(ylims)
ax.text(-0.15, 1.15, panelnames[6], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("p=0.1")
xlabel("τ")
grid()
