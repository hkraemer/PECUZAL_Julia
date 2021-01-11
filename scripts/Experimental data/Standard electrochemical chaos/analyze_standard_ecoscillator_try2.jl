using DrWatson
@quickactivate "PECUZAL_Julia"

#using DelayEmbeddings
using DelimitedFiles
using StatsBase
using Plots
using StatsPlots
using Statistics

include("../../../src/pecuzal_method.jl")
include("../../../src/data_analysis_functions.jl")

## We analyze the reconstruction from standard time delay embedding, the
# MDOP embedding algorithm and the Garcia & Almeida method and compare it to
# our propsed method: The pecuzal method.
# For evaluation of the different reconstructions we consider Uzal's L-statistic,
# the fraction of recurrence rates from the JRP of the reconstruction and the
# reference and only the reference. As a third criterion we look at the
# recurrence time entropy, since it is related to the Kolmogorov-entropy as a
# dynamical invariant. We also look at two other RQA-measures and, moreover, the
# generalized mutual false nearest neighbors.

# Here we look at real data of a chaotic electrochemical oscillator (Ni) from
# the Kiss group with two different resistances R₁ = 1,5kΩ, R₂ = 2,5kΩ leading to
# "more" chaoticity.

dim_tde_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_tde_ref1.csv")
dim_tde_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_tde_ref2.csv")
dim_GA_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_GA_ref1.csv")
dim_GA_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_GA_ref2.csv")
dim_mdop_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_mdop_ref1.csv")
dim_mdop_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_mdop_ref2.csv")
dim_pec_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_pec_ref1.csv")
dim_pec_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_pec_ref2.csv")

Tw_tde_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_tde_ref1.csv")
Tw_tde_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_tde_ref2.csv")
Tw_GA_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_GA_ref1.csv")
Tw_GA_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_GA_ref2.csv")
Tw_mdop_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_mdop_ref1.csv")
Tw_mdop_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_mdop_ref2.csv")
Tw_pec_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_pec_ref1.csv")
Tw_pec_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_pec_ref2.csv")

RQA_tde_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_tde_ref1.csv")
RQA_tde_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_tde_ref2.csv")
RQA_GA_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_GA_ref1.csv")
RQA_GA_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_GA_ref2.csv")
RQA_mdop_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_mdop_ref1.csv")
RQA_mdop_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_mdop_ref2.csv")
RQA_pec_ref1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_pec_ref1.csv")
RQA_pec_ref2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_pec_ref2.csv")

RQA_tde1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_tde1.csv")
RQA_tde2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_tde2.csv")
RQA_GA1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_GA1.csv")
RQA_GA2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_GA2.csv")
RQA_mdop1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_mdop1.csv")
RQA_mdop2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_mdop2.csv")
RQA_pec1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_pec1.csv")
RQA_pec2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_pec2.csv")

dim_tde1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_tde1.csv")
dim_tde2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_tde2.csv")
dim_GA1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_GA1.csv")
dim_GA2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_GA2.csv")
dim_mdop1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_mdop1.csv")
dim_mdop2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_mdop2.csv")
dim_pec1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_pec1.csv")
dim_pec2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_pec2.csv")

Tw_tde1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_tde1.csv")
Tw_tde2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_tde2.csv")
Tw_GA1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_GA1.csv")
Tw_GA2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_GA2.csv")
Tw_mdop1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_mdop1.csv")
Tw_mdop2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_mdop2.csv")
Tw_pec1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_pec1.csv")
Tw_pec2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/Tw_pec2.csv")

# basic statistics
mean_RQA_tde1 = mean(RQA_tde1; dims=2)
mean_RQA_tde2 = mean(RQA_tde2; dims=2)
median_RQA_tde1 = median(RQA_tde1; dims=2)
median_RQA_tde2 = median(RQA_tde2; dims=2)
std_RQA_tde1 = std(RQA_tde1; dims=2)
std_RQA_tde2 = std(RQA_tde2; dims=2)
mean_Tw_tde1 = mean(Tw_tde1)
mean_Tw_tde2 = mean(Tw_tde2)
std_Tw_tde1 = std(Tw_tde1)
std_Tw_tde2 = std(Tw_tde2)
mean_dim_tde1 = mean(dim_tde1)
mean_dim_tde2 = mean(dim_tde2)
std_dim_tde1 = std(dim_tde1)
std_dim_tde2 = std(dim_tde2)

mean_RQA_GA1 = mean(RQA_GA1; dims=2)
mean_RQA_GA2 = mean(RQA_GA2; dims=2)
median_RQA_GA1 = median(RQA_GA1; dims=2)
median_RQA_GA2 = median(RQA_GA2; dims=2)
std_RQA_GA1 = std(RQA_GA1; dims=2)
std_RQA_GA2 = std(RQA_GA2; dims=2)
mean_Tw_GA1 = mean(Tw_GA1)
mean_Tw_GA2 = mean(Tw_GA2)
std_Tw_GA1 = std(Tw_GA1)
std_Tw_GA2 = std(Tw_GA2)
mean_dim_GA1 = mean(dim_GA1)
mean_dim_GA2 = mean(dim_GA2)
std_dim_GA1 = std(dim_GA1)
std_dim_GA2 = std(dim_GA2)

mean_RQA_mdop1 = mean(RQA_mdop1; dims=2)
mean_RQA_mdop2 = mean(RQA_mdop2; dims=2)
median_RQA_mdop1 = median(RQA_mdop1; dims=2)
median_RQA_mdop2 = median(RQA_mdop2; dims=2)
std_RQA_mdop1 = std(RQA_mdop1; dims=2)
std_RQA_mdop2 = std(RQA_mdop2; dims=2)
mean_Tw_mdop1 = mean(Tw_mdop1)
mean_Tw_mdop2 = mean(Tw_mdop2)
std_Tw_mdop1 = std(Tw_mdop1)
std_Tw_mdop2 = std(Tw_mdop2)
mean_dim_mdop1 = mean(dim_mdop1)
mean_dim_mdop2 = mean(dim_mdop2)
std_dim_mdop1 = std(dim_mdop1)
std_dim_mdop2 = std(dim_mdop2)

mean_RQA_pec1 = mean(RQA_pec1; dims=2)
mean_RQA_pec2 = mean(RQA_pec2; dims=2)
median_RQA_pec1 = median(RQA_pec1; dims=2)
median_RQA_pec2 = median(RQA_pec2; dims=2)
std_RQA_pec1 = std(RQA_pec1; dims=2)
std_RQA_pec2 = std(RQA_pec2; dims=2)
mean_Tw_pec1 = mean(Tw_pec1)
mean_Tw_pec2 = mean(Tw_pec2)
std_Tw_pec1 = std(Tw_pec1)
std_Tw_pec2 = std(Tw_pec2)
mean_dim_pec1 = mean(dim_pec1)
mean_dim_pec2 = mean(dim_pec2)
std_dim_pec1 = std(dim_pec1)
std_dim_pec2 = std(dim_pec2)


# compute histograms

p_tde_ENTR1 = histogram(RQA_tde1[1,:],title = "TDE 1 ENTR")
p_tde_ENTR2 = histogram(RQA_tde2[1,:],title = "TDE 2 ENTR")
p_tde_LAM1 = histogram(RQA_tde1[2,:],title = "TDE 1 LAM")
p_tde_LAM2 = histogram(RQA_tde2[2,:],title = "TDE 2 LAM")
p_tde_RTE1 = histogram(RQA_tde1[3,:],title = "TDE 1 RTE")
p_tde_RTE2 = histogram(RQA_tde2[3,:],title = "TDE 2 RTE")
p_tde_TRANS1 = histogram(RQA_tde1[4,:],title = "TDE 1 TRANS")
p_tde_TRANS2 = histogram(RQA_tde2[4,:],title = "TDE 2 TRANS")

p_GA_ENTR1 = histogram(RQA_GA1[1,:],title = "GA 1 ENTR")
p_GA_ENTR2 = histogram(RQA_GA2[1,:],title = "GA 2 ENTR")
p_GA_LAM1 = histogram(RQA_GA1[2,:],title = "GA 1 LAM")
p_GA_LAM2 = histogram(RQA_GA2[2,:],title = "GA 2 LAM")
p_GA_RTE1 = histogram(RQA_GA1[3,:],title = "GA 1 RTE")
p_GA_RTE2 = histogram(RQA_GA2[3,:],title = "GA 2 RTE")
p_GA_TRANS1 = histogram(RQA_GA1[4,:],title = "GA 1 TRANS")
p_GA_TRANS2 = histogram(RQA_GA2[4,:],title = "GA 2 TRANS")

p_mdop_ENTR1 = histogram(RQA_mdop1[1,:],title = "mdop 1 ENTR")
p_mdop_ENTR2 = histogram(RQA_mdop2[1,:],title = "mdop 2 ENTR")
p_mdop_LAM1 = histogram(RQA_mdop1[2,:],title = "mdop 1 LAM")
p_mdop_LAM2 = histogram(RQA_mdop2[2,:],title = "mdop 2 LAM")
p_mdop_RTE1 = histogram(RQA_mdop1[3,:],title = "mdop 1 RTE")
p_mdop_RTE2 = histogram(RQA_mdop2[3,:],title = "mdop 2 RTE")
p_mdop_TRANS1 = histogram(RQA_mdop1[4,:],title = "mdop 1 TRANS")
p_mdop_TRANS2 = histogram(RQA_mdop2[4,:],title = "mdop 2 TRANS")

p_pec_ENTR1 = histogram(RQA_pec1[1,:],title = "pec 1 ENTR")
p_pec_ENTR2 = histogram(RQA_pec2[1,:],title = "pec 2 ENTR")
p_pec_LAM1 = histogram(RQA_pec1[2,:],title = "pec 1 LAM")
p_pec_LAM2 = histogram(RQA_pec2[2,:],title = "pec 2 LAM")
p_pec_RTE1 = histogram(RQA_pec1[3,:],title = "pec 1 RTE")
p_pec_RTE2 = histogram(RQA_pec2[3,:],title = "pec 2 RTE")
p_pec_TRANS1 = histogram(RQA_pec1[4,:],title = "pec 1 TRANS")
p_pec_TRANS2 = histogram(RQA_pec2[4,:],title = "pec 2 TRANS")


## ENTR

plot(p_tde_ENTR1, label="data from subsamples")
plot!([RQA_tde_ref1[1]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_tde1[1]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_tde_ENTR2, label="data from subsamples")
plot!([RQA_tde_ref2[1]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_tde2[1]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_GA_ENTR1, label="data from subsamples")
plot!([RQA_GA_ref1[1]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_GA1[1]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_GA_ENTR2, label="data from subsamples")
plot!([RQA_GA_ref2[1]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_GA2[1]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_mdop_ENTR1, label="data from subsamples")
plot!([RQA_mdop_ref1[1]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_mdop1[1]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_mdop_ENTR2, label="data from subsamples")
plot!([RQA_mdop_ref2[1]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_mdop2[1]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_pec_ENTR1, label="data from subsamples")
plot!([RQA_pec_ref1[1]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_pec1[1]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_pec_ENTR2, label="data from subsamples")
plot!([RQA_pec_ref2[1]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_pec2[1]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")


## LAM

plot(p_tde_LAM1, label="data from subsamples")
plot!([RQA_tde_ref1[2]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_tde1[2]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_tde_LAM2, label="data from subsamples")
plot!([RQA_tde_ref2[2]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_tde2[2]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_GA_LAM1, label="data from subsamples")
plot!([RQA_GA_ref1[2]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_GA1[2]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_GA_LAM2, label="data from subsamples")
plot!([RQA_GA_ref2[2]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_GA2[2]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_mdop_LAM1, label="data from subsamples")
plot!([RQA_mdop_ref1[2]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_mdop1[2]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_mdop_LAM2, label="data from subsamples")
plot!([RQA_mdop_ref2[2]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_mdop2[2]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_pec_LAM1, label="data from subsamples")
plot!([RQA_pec_ref1[2]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_pec1[2]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_pec_LAM2, label="data from subsamples")
plot!([RQA_pec_ref2[2]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_pec2[2]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

## RTE

plot(p_tde_RTE1, label="data from subsamples")
plot!([RQA_tde_ref1[3]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_tde1[3]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_tde_RTE2, label="data from subsamples")
plot!([RQA_tde_ref2[3]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_tde2[3]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_GA_RTE1, label="data from subsamples")
plot!([RQA_GA_ref1[3]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_GA1[3]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_GA_RTE2, label="data from subsamples")
plot!([RQA_GA_ref2[3]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_GA2[3]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_mdop_RTE1, label="data from subsamples")
plot!([RQA_mdop_ref1[3]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_mdop1[3]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_mdop_RTE2, label="data from subsamples")
plot!([RQA_mdop_ref2[3]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_mdop2[3]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_pec_RTE1, label="data from subsamples")
plot!([RQA_pec_ref1[3]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_pec1[3]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_pec_RTE2, label="data from subsamples")
plot!([RQA_pec_ref2[3]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_pec2[3]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

## TRANS

plot(p_tde_TRANS1, label="data from subsamples")
plot!([RQA_tde_ref1[4]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_tde1[4]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_tde_TRANS2, label="data from subsamples")
plot!([RQA_tde_ref2[4]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_tde2[4]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_GA_TRANS1, label="data from subsamples")
plot!([RQA_GA_ref1[4]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_GA1[4]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_GA_TRANS2, label="data from subsamples")
plot!([RQA_GA_ref2[4]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_GA2[4]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_mdop_TRANS1, label="data from subsamples")
plot!([RQA_mdop_ref1[4]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_mdop1[4]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_mdop_TRANS2, label="data from subsamples")
plot!([RQA_mdop_ref2[4]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_mdop2[4]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_pec_TRANS1, label="data from subsamples")
plot!([RQA_pec_ref1[4]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_pec1[4]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

plot(p_pec_TRANS2, label="data from subsamples")
plot!([RQA_pec_ref2[4]], seriestype=vline, linecolor="red", linewidth=10, label="reference")
plot!([median_RQA_pec2[4]], seriestype=vline, linecolor="gray", linestyle= :dash, linewidth=10, label="median")

## Deviation of median from reference in σ-units:

diff_RQA_tde1 = zeros(4)
diff_RQA_tde2 = zeros(4)
diff_RQA_GA1 = zeros(4)
diff_RQA_GA2 = zeros(4)
diff_RQA_mdop1 = zeros(4)
diff_RQA_mdop2 = zeros(4)
diff_RQA_pec1 = zeros(4)
diff_RQA_pec2 = zeros(4)

for i = 1:4
    diff_RQA_tde1[i] = abs(median_RQA_tde1[i]-RQA_tde_ref1[i])/std_RQA_tde1[i]
    diff_RQA_tde2[i] = abs(median_RQA_tde2[i]-RQA_tde_ref2[i])/std_RQA_tde2[i]

    diff_RQA_GA1[i] = abs(median_RQA_GA1[i]-RQA_GA_ref1[i])/std_RQA_GA1[i]
    diff_RQA_GA2[i] = abs(median_RQA_GA2[i]-RQA_GA_ref2[i])/std_RQA_GA2[i]

    diff_RQA_mdop1[i] = abs(median_RQA_mdop1[i]-RQA_mdop_ref1[i])/std_RQA_mdop1[i]
    diff_RQA_mdop2[i] = abs(median_RQA_mdop2[i]-RQA_mdop_ref2[i])/std_RQA_mdop2[i]

    diff_RQA_pec1[i] = abs(median_RQA_pec1[i]-RQA_pec_ref1[i])/std_RQA_pec1[i]
    diff_RQA_pec2[i] = abs(median_RQA_pec2[i]-RQA_pec_ref2[i])/std_RQA_pec2[i]
end

rel_diff_RQA_tde1 = zeros(4)
rel_diff_RQA_tde2 = zeros(4)
rel_diff_RQA_GA1 = zeros(4)
rel_diff_RQA_GA2 = zeros(4)
rel_diff_RQA_mdop1 = zeros(4)
rel_diff_RQA_mdop2 = zeros(4)
rel_diff_RQA_pec1 = zeros(4)
rel_diff_RQA_pec2 = zeros(4)

for i = 1:4
    rel_diff_RQA_tde1[i] = abs(median_RQA_tde1[i]-RQA_tde_ref1[i])/RQA_tde_ref1[i]
    rel_diff_RQA_tde2[i] = abs(median_RQA_tde2[i]-RQA_tde_ref2[i])/RQA_tde_ref2[i]

    rel_diff_RQA_GA1[i] = abs(median_RQA_GA1[i]-RQA_GA_ref1[i])/RQA_GA_ref1[i]
    rel_diff_RQA_GA2[i] = abs(median_RQA_GA2[i]-RQA_GA_ref2[i])/RQA_GA_ref2[i]

    rel_diff_RQA_mdop1[i] = abs(median_RQA_mdop1[i]-RQA_mdop_ref1[i])/RQA_mdop_ref1[i]
    rel_diff_RQA_mdop2[i] = abs(median_RQA_mdop2[i]-RQA_mdop_ref2[i])/RQA_mdop_ref2[i]

    rel_diff_RQA_pec1[i] = abs(median_RQA_pec1[i]-RQA_pec_ref1[i])/RQA_pec_ref1[i]
    rel_diff_RQA_pec2[i] = abs(median_RQA_pec2[i]-RQA_pec_ref2[i])/RQA_pec_ref2[i]
end


mn1 = [diff_RQA_tde1[1], diff_RQA_tde1[2], diff_RQA_tde1[3], diff_RQA_tde1[4],
    diff_RQA_GA1[1], diff_RQA_GA1[2], diff_RQA_GA1[3], diff_RQA_GA1[4],
    diff_RQA_mdop1[1], diff_RQA_mdop1[2], diff_RQA_mdop1[3], diff_RQA_mdop1[4],
    diff_RQA_pec1[1], diff_RQA_pec1[2], diff_RQA_pec1[3], diff_RQA_pec1[4]
    ]
mn2 = [diff_RQA_tde2[1], diff_RQA_tde2[2], diff_RQA_tde2[3], diff_RQA_tde2[4],
    diff_RQA_GA2[1], diff_RQA_GA2[2], diff_RQA_GA2[3], diff_RQA_GA2[4],
    diff_RQA_mdop2[1], diff_RQA_mdop2[2], diff_RQA_mdop2[3], diff_RQA_mdop2[4],
    diff_RQA_pec2[1], diff_RQA_pec2[2], diff_RQA_pec2[3], diff_RQA_pec2[4]
    ]
sx = repeat(["Cao's TDE", "Garcia & Almeida", "MDOP", "PECUZAL"], inner = 4)
nam = repeat(["ENTR", "LAM", "RTE", "TRANS"], outer=4)


StatsPlots.groupedbar(nam, mn1, group = sx, ylabel = "deviation from reference [σ]",
        title = "chem. Oscillators Exp.1")

StatsPlots.groupedbar(nam, mn2, group = sx, ylabel = "deviation from reference [σ]",
        title = "chem. Oscillators Exp.2")

total_diff_tde1 = sum(diff_RQA_tde1)
total_diff_tde2 = sum(diff_RQA_tde2)

total_diff_GA1 = sum(diff_RQA_GA1)
total_diff_GA2 = sum(diff_RQA_GA2)

total_diff_mdop1 = sum(diff_RQA_mdop1)
total_diff_mdop2 = sum(diff_RQA_mdop2)

total_diff_pec1 = sum(diff_RQA_pec1)
total_diff_pec2 = sum(diff_RQA_pec2)


rel_mn1 = [rel_diff_RQA_tde1[1], rel_diff_RQA_tde1[2], rel_diff_RQA_tde1[3], rel_diff_RQA_tde1[4],
    rel_diff_RQA_GA1[1], rel_diff_RQA_GA1[2], rel_diff_RQA_GA1[3], rel_diff_RQA_GA1[4],
    rel_diff_RQA_mdop1[1], rel_diff_RQA_mdop1[2], rel_diff_RQA_mdop1[3], rel_diff_RQA_mdop1[4],
    rel_diff_RQA_pec1[1], rel_diff_RQA_pec1[2], rel_diff_RQA_pec1[3], rel_diff_RQA_pec1[4]
    ]
rel_mn2 = [rel_diff_RQA_tde2[1], rel_diff_RQA_tde2[2], rel_diff_RQA_tde2[3], rel_diff_RQA_tde2[4],
    rel_diff_RQA_GA2[1], rel_diff_RQA_GA2[2], rel_diff_RQA_GA2[3], rel_diff_RQA_GA2[4],
    rel_diff_RQA_mdop2[1], rel_diff_RQA_mdop2[2], rel_diff_RQA_mdop2[3], rel_diff_RQA_mdop2[4],
    rel_diff_RQA_pec2[1], rel_diff_RQA_pec2[2], rel_diff_RQA_pec2[3], rel_diff_RQA_pec2[4]
    ]

# save relative deviations for plotting in Matlab
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/relative_dev_RQA_1.csv", rel_mn1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/relative_dev_RQA_2.csv", rel_mn2)

StatsPlots.groupedbar(nam, rel_mn1, group = sx, ylabel = "rel .deviation from reference",
        title = "chem. Oscillators Exp.1")

StatsPlots.groupedbar(nam, rel_mn2, group = sx, ylabel = "rel .deviation from reference",
        title = "chem. Oscillators Exp.2")

total_rel_diff_tde1 = sum(rel_diff_RQA_tde1)
total_rel_diff_tde2 = sum(rel_diff_RQA_tde2)

total_rel_diff_GA1 = sum(rel_diff_RQA_GA1)
total_rel_diff_GA2 = sum(rel_diff_RQA_GA2)

total_rel_diff_mdop1 = sum(rel_diff_RQA_mdop1)
total_rel_diff_mdop2 = sum(rel_diff_RQA_mdop2)

total_rel_diff_pec1 = sum(rel_diff_RQA_pec1)
total_rel_diff_pec2 = sum(rel_diff_RQA_pec2)

# for plotting 1-rel.diff:
plot_rel_mn1 = [1-rel_diff_RQA_tde1[1], 1-rel_diff_RQA_tde1[2], 1-rel_diff_RQA_tde1[3], 1-rel_diff_RQA_tde1[4],
    1-rel_diff_RQA_GA1[1], 1-rel_diff_RQA_GA1[2], 1-rel_diff_RQA_GA1[3], 1-rel_diff_RQA_GA1[4],
    1-rel_diff_RQA_mdop1[1], 1-rel_diff_RQA_mdop1[2], 1-rel_diff_RQA_mdop1[3], 1-rel_diff_RQA_mdop1[4],
    1-rel_diff_RQA_pec1[1], 1-rel_diff_RQA_pec1[2], 1-rel_diff_RQA_pec1[3], 1-rel_diff_RQA_pec1[4]
    ]
plot_rel_mn2 = [1-rel_diff_RQA_tde2[1], 1-rel_diff_RQA_tde2[2], 1-rel_diff_RQA_tde2[3], 1-rel_diff_RQA_tde2[4],
    1-rel_diff_RQA_GA2[1], 1-rel_diff_RQA_GA2[2], 1-rel_diff_RQA_GA2[3], 1-rel_diff_RQA_GA2[4],
    1-rel_diff_RQA_mdop2[1], 1-rel_diff_RQA_mdop2[2], 1-rel_diff_RQA_mdop2[3], 1-rel_diff_RQA_mdop2[4],
    1-rel_diff_RQA_pec2[1], 1-rel_diff_RQA_pec2[2], 1-rel_diff_RQA_pec2[3], 1-rel_diff_RQA_pec2[4]
    ]

StatsPlots.groupedbar(nam, plot_rel_mn1, group = sx, ylabel = "rel .deviation from reference",
        title = "chem. Oscillators Exp.1")

StatsPlots.groupedbar(nam, plot_rel_mn2, group = sx, ylabel = "rel .deviation from reference",
        title = "chem. Oscillators Exp.2")


# display results:
round.(plot_rel_mn1[5]; digits=4)
