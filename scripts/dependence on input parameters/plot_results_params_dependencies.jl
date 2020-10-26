using DrWatson
@quickactivate "new-embedding-methods"

using DelimitedFiles, PyPlot
pygui(true)

# Here plot the results from scripts `dependence_on_alpha.jl`, `dependence_on_K.jl`,
# `dependence_on_p.jl`, `dependence_on_delta_neighborhood.jl`, `dependence_on_Tw.jl`,

panelnames = ["A" "B" "C" "D" "E" "F" "G" "H" "I" "J"]

fig = figure(figsize=[15,12])

lwg = 3         # linewidth of the graph
lwa = 2         # linewidth of the axis
fsp = 20        # Fontsize of panelnames
fsa = 14        # Fontsize of the axis
fsl = 12        # Fontsize of the legendentries
fst = 16        # Fontsize of title
axislabelsize = 14 # axislabelsize
ticklabelsize = 12  # labelsize of ticks
ms = 10         # markersize of maxima
ms_style = "*"  # markerstyle

subplots_adjust(left = 0.07,
right = 0.96,
bottom = 0.05,
top = 0.96,
wspace = 0.25,
hspace = 0.4)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
font0 = Dict(
        "font.size" => fsa,
        "font.weight" => "bold",
        "axes.labelweight" => "bold",
        "axes.titleweight" => "bold",
        "axes.labelsize" => axislabelsize,
        "axes.linewidth" => lwa,
        "xtick.labelsize" => ticklabelsize,
        "ytick.labelsize" => ticklabelsize,
        "legend.fontsize" => fsl)
merge!(rcParams, font0)

# 1 Dependence on KNN

Ks = readdlm("./scripts/computed data/dependence_on_K_Ks.csv")
sizes = readdlm("./scripts/computed data/dependence_on_K_sizes.csv")
Ls = readdlm("./scripts/computed data/dependence_on_K_Ls.csv")
τ_vals = readdlm("./scripts/computed data/dependence_on_K_tau_vals.csv")
τ_vals = Int.(τ_vals)

ax = fig.add_subplot(5, 2, 1)
p1 = plot(Ks, Ls[:,1], label="1st embedding cylce", linewidth=lwg)
color1 = p1[1].get_color()
scatter(Ks, Ls[:,1], label="")
p2 = plot(Ks, Ls[:,2], label="2nd embedding cylce", linewidth=lwg)
color2 = p2[1].get_color()
scatter(Ks, Ls[:,2], label="")
legend()
ax.text(-0.15, 1.1, panelnames[1], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("L statistics")
xlabel("K-NN")
ylabel("L")
grid()

ax = fig.add_subplot(5, 2, 2)
plot(Ks, τ_vals[:,1], label="τ₀, 1st embedding cylce", linewidth=lwg, color="r")
scatter(Ks, τ_vals[:,1], label="", color="r")
plot(Ks, τ_vals[:,2], label="τ₁, 2nd embedding cylce", linewidth=lwg, color=color1)
scatter(Ks, τ_vals[:,2], label="", color=color1)
plot(Ks, τ_vals[:,3], label="τ₂, 3rd embedding cylce", linewidth=lwg, color=color2)
scatter(Ks, τ_vals[:,3], label="", color=color2)
legend(loc=1)
ax.set_ylim([-1,24])
ax.text(-0.15, 1.1, panelnames[2], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("chosen delays")
xlabel("k-NN")
ylabel("τ")
grid()

# 2) Dependence on δ-Neighborhood size

deltas = readdlm("./scripts/computed data/dependence_on_delta_neighborhood_Ks.csv")
sizes = readdlm("./scripts/computed data/dependence_on_delta_neighborhood_sizes.csv")
Ls = readdlm("./scripts/computed data/dependence_on_delta_neighborhood_Ls.csv")
τ_vals = readdlm("./scripts/computed data/dependence_on_delta_neighborhood_tau_vals.csv")
τ_vals = Int.(τ_vals)

ax = fig.add_subplot(5, 2, 3)
p1 = plot(deltas, Ls[:,1], label="1st embedding cylce", linewidth=lwg)
color1 = p1[1].get_color()
scatter(deltas, Ls[:,1], label="")
p2 = plot(deltas, Ls[:,2], label="2nd embedding cylce", linewidth=lwg)
color2 = p2[1].get_color()
scatter(deltas, Ls[:,2], label="")
#legend()
ax.text(-0.15, 1.1, panelnames[3], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
#title("minimum L")
xlabel(L"\delta-Neighborhood size", fontweight="normal")
ylabel("L")
grid()

ax = fig.add_subplot(5, 2, 4)
plot(deltas, τ_vals[:,1], label="τ₁, 1st embedding cylce", linewidth=lwg, color="r")
scatter(deltas, τ_vals[:,1], label="", color="r")
plot(deltas, τ_vals[:,2], label="τ₂, 2nd embedding cylce", linewidth=lwg, color=color1)
scatter(deltas, τ_vals[:,2], label="", color=color1)
plot(deltas, τ_vals[:,3], label="τ₃, 3rd embedding cylce", linewidth=lwg, color=color2)
scatter(deltas, τ_vals[:,3], label="", color=color2)
#legend()
ax.set_ylim([-1,24])
ax.text(-0.15, 1.1, panelnames[4], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
#title("chosen delays")
xlabel(L"\delta-Neighborhood size")
ylabel("τ")
grid()

# 3) Dependence on p

ps = readdlm("./scripts/computed data/dependence_on_p_ps.csv")
sizes = readdlm("./scripts/computed data/dependence_on_p_sizes.csv")
Ls = readdlm("./scripts/computed data/dependence_on_p_Ls.csv")
τ_vals = readdlm("./scripts/computed data/dependence_on_p_tau_vals.csv")
τ_vals = Int.(τ_vals)
εs = readdlm("./scripts/computed data/dependence_on_p_epsilons.csv")


ax = fig.add_subplot(5, 2, 5)
p1 = plot(ps', Ls[:,1], label="1st embedding cylce", linewidth=lwg)
color1 = p1[1].get_color()
scatter(ps', Ls[:,1], label="")
p2 = plot(ps', Ls[:,2], label="2nd embedding cylce", linewidth=lwg)
color2 = p2[1].get_color()
scatter(ps', Ls[:,2], label="")
#legend()
ax.text(-0.15, 1.1, panelnames[5], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
#title("minimum L")
xlabel("(binomial-) p")
ylabel("L")
grid()

ax = fig.add_subplot(5, 2, 6)
plot(ps', τ_vals[:,1], label="τ₁, 1st embedding cylce", linewidth=lwg, color="r")
scatter(ps', τ_vals[:,1], label="", color="r")
plot(ps', τ_vals[:,2], label="τ₂, 2nd embedding cylce", linewidth=lwg, color=color1)
scatter(ps', τ_vals[:,2], label="", color=color1)
plot(ps', τ_vals[:,3], label="τ₃, 3rd embedding cylce", linewidth=lwg, color=color2)
scatter(ps', τ_vals[:,3], label="", color=color2)
#legend()
ax.set_ylim([-1,24])
ax.text(-0.15, 1.1, panelnames[6], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
#title("chosen delays")
xlabel("(binomial-) p")
ylabel("τ")
grid()


# 4) Dependence on α

alphas = readdlm("./scripts/computed data/dependence_on_alpha_alphas.csv")
sizes = readdlm("./scripts/computed data/dependence_on_alpha_sizes.csv")
Ls = readdlm("./scripts/computed data/dependence_on_alpha_Ls.csv")
τ_vals = readdlm("./scripts/computed data/dependence_on_alpha_tau_vals.csv")
τ_vals = Int.(τ_vals)
εs = readdlm("./scripts/computed data/dependence_on_alpha_epsilons.csv")


ax = fig.add_subplot(5, 2, 7)
p1 = plot(alphas', Ls[:,1], label="1st embedding cylce", linewidth=lwg)
color1 = p1[1].get_color()
scatter(alphas', Ls[:,1], label="")
p2 = plot(alphas', Ls[:,2], label="2nd embedding cylce", linewidth=lwg)
color2 = p2[1].get_color()
scatter(alphas', Ls[:,2], label="")
#legend()
ax.text(-0.15, 1.1, panelnames[7], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
ax.set_xticks([0.01, 0.05])
#title("minimum L")
xlabel("(binomial-) α")
ylabel("L")
grid()

ax = fig.add_subplot(5, 2, 8)
plot(alphas', τ_vals[:,1], label="τ₁, embedding cylce 0", linewidth=lwg, color="r")
scatter(alphas', τ_vals[:,1], label="", color="r")
plot(alphas', τ_vals[:,2], label="τ₂, embedding cylce 1", linewidth=lwg, color=color1)
scatter(alphas', τ_vals[:,2], label="", color=color1)
plot(alphas', τ_vals[:,3], label="τ₃, embedding cylce 2", linewidth=lwg, color=color2)
scatter(alphas', τ_vals[:,3], label="", color=color2)
#legend()
ax.set_ylim([-1,24])
ax.set_xticks([0.01, 0.05])
ax.text(-0.15, 1.1, panelnames[8], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
#title("chosen delays")
xlabel("(binomial-) α")
ylabel("τ")
grid()

# 5) Dependence on Tw

Tws = readdlm("./scripts/computed data/dependence_on_Tw_Tws.csv")
sizes = readdlm("./scripts/computed data/dependence_on_Tw_sizes.csv")
Ls = readdlm("./scripts/computed data/dependence_on_Tw_Ls.csv")
τ_vals = readdlm("./scripts/computed data/dependence_on_Tw_tau_vals.csv")

τ_vals[1,3]=NaN64
τ_vals[2,3]=NaN64
τ_vals[3,3]=NaN64
Ls[1,3]=NaN64
Ls[2,3]=NaN64
Ls[3,3]=NaN64

ax = fig.add_subplot(5, 2, 9)
p1 = plot(Tws, Ls[:,1], label="embedding cylce 1", linewidth=lwg)
color1 = p1[1].get_color()
scatter(Tws, Ls[:,1], label="")
p2 = plot(Tws, Ls[:,2], label="embedding cylce 2", linewidth=lwg)
color2 = p2[1].get_color()
scatter(Tws, Ls[:,2], label="")
#legend()
ax.text(-0.15, 1.1, panelnames[9], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
#title("minimum L")
xlabel("time window")
ylabel("L")
grid()

ax = fig.add_subplot(5, 2, 10)
plot(Tws, τ_vals[:,1], label="τ₁, embedding cylce 0", linewidth=lwg, color="r")
scatter(Tws, τ_vals[:,1], label="", color="r")
plot(Tws, τ_vals[:,2], label="τ₂, embedding cylce 1", linewidth=lwg, color=color1)
scatter(Tws, τ_vals[:,2], label="", color=color1)
plot(Tws, τ_vals[:,3], label="τ₃, embedding cylce 2", linewidth=lwg, color=color2)
scatter(Tws, τ_vals[:,3], label="", color=color2)
#legend()
ax.set_ylim([-1,24])
ax.text(-0.15, 1.1, panelnames[10], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
#title("chosen delays")
xlabel("time window")
ylabel("τ")
grid()




## Figure for ε*-statistics with varying p

ps = readdlm("./scripts/computed data/dependence_on_p_ps.csv")
sizes = readdlm("./scripts/computed data/dependence_on_p_sizes.csv")
Ls = readdlm("./scripts/computed data/dependence_on_p_Ls.csv")
τ_vals = readdlm("./scripts/computed data/dependence_on_p_tau_vals.csv")
τ_vals = Int.(τ_vals)
εs = readdlm("./scripts/computed data/dependence_on_p_epsilons.csv")

Tmax=100
taus = vec(0:Tmax)
NN =length(taus)

εs_ = [zeros(NN,3), zeros(NN,3), zeros(NN,3), zeros(NN,3), zeros(NN,3), zeros(NN,3)]
for j = 1:6
        for i = 1:3
                εs_[j][:,i] = εs[j,(i-1)*100+i:i*NN]
        end
end
εs=εs_

panelnames = ["A" "B" "C" "D" "E" "F"]
titlenames = ["p=0.6" "p=0.5" "p=0.4" "p=0.3" "p=0.2" "p=0.1"]

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

for i = 1:6
        ax = fig.add_subplot(2, 3, i)
        plot(taus,εs[i][:,1], label="τ₀=$(τ_vals[i,1])", linewidth=lwg)
        plot(taus,εs[i][:,2], label="τ₁=$(τ_vals[i,2])", linewidth=lwg)
        plot(taus,εs[i][:,3], label="τ₂=$(τ_vals[i,3])", linewidth=lwg)
        plot([taus[τ_vals[i,2]+1]], [εs[i][τ_vals[i,2]+1,1]], ms_style, label="", markersize = ms, color = "black")
        plot([taus[τ_vals[i,3]+1]], [εs[i][τ_vals[i,3]+1,2]], ms_style, label="", markersize = ms, color = "black")
        legend()
        ax.set_ylim(ylims)
        ax.text(-0.15, 1.15, panelnames[i], transform=ax.transAxes,
        fontsize=fsp, fontweight="bold", va="top")
        title(titlenames[i])
        if i > 3
        xlabel("τ")
        end
        if i == 1
        ylabel("⟨ε★⟩")
        end
        grid()
end


## Figure for ε*-statistics with varying α

alphas = readdlm("./scripts/computed data/dependence_on_alpha_alphas.csv")
sizes = readdlm("./scripts/computed data/dependence_on_alpha_sizes.csv")
Ls = readdlm("./scripts/computed data/dependence_on_alpha_Ls.csv")
τ_vals = readdlm("./scripts/computed data/dependence_on_alpha_tau_vals.csv")
τ_vals = Int.(τ_vals)
εs = readdlm("./scripts/computed data/dependence_on_alpha_epsilons.csv")


Tmax=100
taus = vec(0:Tmax)
NN =length(taus)

εs_ = [zeros(NN,3), zeros(NN,3)]
for j = 1:2
        for i = 1:3
                εs_[j][:,i] = εs[j,(i-1)*100+i:i*NN]
        end
end
εs=εs_

τ_vals = Int.(τ_vals)

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
plot(taus,εs[1][:,1], label="τ₀=$(τ_vals[1,1])", linewidth = lwg)
plot(taus,εs[1][:,2], label="τ₁=$(τ_vals[1,2])", linewidth = lwg)
plot(taus,εs[1][:,3], label="τ₂=$(τ_vals[1,3])", linewidth = lwg)
plot([taus[τ_vals[1,2]+1]], [εs[1][τ_vals[1,2]+1,1]], ms_style, label="", markersize = ms, color = "black")
plot([taus[τ_vals[1,3]+1]], [εs[1][τ_vals[1,3]+1,2]], ms_style, label="", markersize = ms, color = "black")
legend()
ax.set_ylim(ylims)
ax.text(-0.15, 1.1, panelnames[1], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("α=0.05")
xlabel("τ")
ylabel("⟨ε★⟩")
grid()

ax = fig.add_subplot(1, 2, 2)
plot(taus,εs[2][:,1], label="τ₀=$(τ_vals[2,1])", linewidth = lwg)
plot(taus,εs[2][:,2], label="τ₁=$(τ_vals[2,2])", linewidth = lwg)
plot(taus,εs[2][:,3], label="τ₂=$(τ_vals[2,3])", linewidth = lwg)
plot([taus[τ_vals[2,2]+1]], [εs[2][τ_vals[2,2]+1,1]], ms_style, label="", markersize = ms, color = "black")
plot([taus[τ_vals[2,3]+1]], [εs[2][τ_vals[2,3]+1,2]], ms_style, label="", markersize = ms, color = "black")
legend()
ax.set_ylim(ylims)
ax.text(-0.15, 1.1, panelnames[2], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("α=0.01")
xlabel("τ")
grid()
