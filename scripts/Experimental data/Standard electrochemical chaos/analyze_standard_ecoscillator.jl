using DrWatson
@quickactivate "PECUZAL_Julia"

using DelayEmbeddings
using DelimitedFiles
using StatsBase
using HypothesisTests

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
L_tde1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_tde1.csv")
L_tde2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_tde2.csv")
dim_GA1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_GA1.csv")
dim_GA2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_GA2.csv")
L_GA1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_GA1.csv")
L_GA2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_GA2.csv")
dim_mdop1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_mdop1.csv")
dim_mdop2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_mdop2.csv")
L_mdop1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_mdop1.csv")
L_mdop2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_mdop2.csv")
dim_pec1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_pec1.csv")
dim_pec2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_pec2.csv")
L_pec1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_pec1.csv")
L_pec2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_pec2.csv")

using Plots

p1 = histogram(RQA_tde1[1,:],title = "TDE 1 ENTR")
p2 = histogram(RQA_GA1[1,:],title = "GA 1 ENTR")
p3 = histogram(RQA_mdop1[1,:],title = "MDOP 1 ENTR")
p4 = histogram(RQA_pec1[1,:],title = "PEC 1 ENTR")
p5 = histogram(RQA_tde2[1,:],title = "TDE 2 ENTR")
p6 = histogram(RQA_GA2[1,:],title = "GA 2 ENTR")
p7 = histogram(RQA_mdop2[1,:],title = "MDOP 2 ENTR")
p8 = histogram(RQA_pec2[1,:],title = "PEC 2 ENTR")
plot(p1,p2,p3,p4,layout = (2, 2), legend = false)
plot(p5,p6,p7,p8,layout = (2, 2), legend = false)

p1 = histogram(RQA_tde1[2,:],title = "TDE 1 TRANS")
p2 = histogram(RQA_GA1[2,:],title = "GA 1 TRANS")
p3 = histogram(RQA_mdop1[2,:],title = "MDOP 1 TRANS")
p4 = histogram(RQA_pec1[2,:],title = "PEC 1 TRANS")
p5 = histogram(RQA_tde2[2,:],title = "TDE 2 TRANS")
p6 = histogram(RQA_GA2[2,:],title = "GA 2 TRANS")
p7 = histogram(RQA_mdop2[2,:],title = "MDOP 2 TRANS")
p8 = histogram(RQA_pec2[2,:],title = "PEC 2 TRANS")
plot(p1,p2,p3,p4,layout = (2, 2), legend = false)
plot(p5,p6,p7,p8,layout = (2, 2), legend = false)

med_TDE = zeros(2,4)
med_GA = zeros(2,4)
med_mdop = zeros(2,4)
med_pec = zeros(2,4)
std_TDE = zeros(2,4)
std_GA = zeros(2,4)
std_mdop = zeros(2,4)
std_pec = zeros(2,4)

Ltde1 = mean(L_tde1)
Lstd_tde1 = std(L_tde1)
Ltde2 = mean(L_tde2)
Lstd_tde2 = std(L_tde2)
dim_tde1 = mean(dim_tde1)
dim_tde2 = mean(dim_tde2)
dimstd_tde1 = std(dim_tde1)
dimstd_tde2 = std(dim_tde2)

LGA1 = mean(L_GA1)
Lstd_GA1 = std(L_GA1)
LGA2 = mean(L_GA2)
Lstd_GA2 = std(L_GA2)
dim_GA1 = mean(dim_GA1)
dim_GA2 = mean(dim_GA2)
dimstd_GA1 = std(dim_GA1)
dimstd_GA2 = std(dim_GA2)

Lmdop1 = mean(L_mdop1)
Lstd_mdop1 = std(L_mdop1)
Lmdop2 = mean(L_mdop2)
Lstd_mdop2 = std(L_mdop2)
dim_mdop1 = mean(dim_mdop1)
dim_mdop2 = mean(dim_mdop2)
dimstd_mdop1 = std(dim_mdop1)
dimstd_mdop2 = std(dim_mdop2)

Lpec1 = mean(L_pec1)
Lstd_pec1 = std(L_pec1)
Lpec2 = mean(L_pec2)
Lstd_pec2 = std(L_pec2)
dim_pec1 = mean(dim_pec1)
dim_pec2 = mean(dim_pec2)
dimstd_pec1 = std(dim_pec1)
dimstd_pec2 = std(dim_pec2)

for i = 1:4
    med_TDE[1,i] = mean(RQA_tde1[i,:])
    med_TDE[2,i] = mean(RQA_tde2[i,:])
    med_GA[1,i] = mean(RQA_GA1[i,:])
    med_GA[2,i] = mean(RQA_GA2[i,:])
    med_mdop[1,i] = mean(RQA_mdop1[i,:])
    med_mdop[2,i] = mean(RQA_mdop2[i,:])
    med_pec[1,i] = mean(RQA_pec1[i,:])
    med_pec[2,i] = mean(RQA_pec2[i,:])

    std_TDE[1,i] = std(RQA_tde1[i,:])
    std_TDE[2,i] = std(RQA_tde2[i,:])
    std_GA[1,i] = std(RQA_GA1[i,:])
    std_GA[2,i] = std(RQA_GA2[i,:])
    std_mdop[1,i] = std(RQA_mdop1[i,:])
    std_mdop[2,i] = std(RQA_mdop2[i,:])
    std_pec[1,i] = std(RQA_pec1[i,:])
    std_pec[2,i] = std(RQA_pec2[i,:])
end

RQA_measure = 4
# 1: ENTR
# 2: LAM
# 3: RTE
# 4: TRANS

wtest_tde = UnequalVarianceTTest(RQA_tde1[RQA_measure,:], RQA_tde2[RQA_measure,:])
wtest_GA = UnequalVarianceTTest(RQA_GA1[RQA_measure,:], RQA_GA2[RQA_measure,:])
wtest_mdop = UnequalVarianceTTest(RQA_mdop1[RQA_measure,:], RQA_mdop2[RQA_measure,:])
wtest_pec = UnequalVarianceTTest(RQA_pec1[RQA_measure,:], RQA_pec2[RQA_measure,:])

pvalue(wtest_tde)
pvalue(wtest_GA)
pvalue(wtest_mdop)
pvalue(wtest_pec)


wtest_tde = UnequalVarianceTTest(L_tde1, L_tde2)
wtest_GA = UnequalVarianceTTest(RQA_GA1[RQA_measure,:], RQA_GA2[RQA_measure,:])
wtest_mdop = UnequalVarianceTTest(RQA_mdop1[RQA_measure,:], RQA_mdop2[RQA_measure,:])
wtest_pec = UnequalVarianceTTest(RQA_pec1[RQA_measure,:], RQA_pec2[RQA_measure,:])
