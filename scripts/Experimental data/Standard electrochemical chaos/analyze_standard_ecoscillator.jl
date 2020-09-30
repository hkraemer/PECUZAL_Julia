using DrWatson
@quickactivate "new-embedding-methods"

using DelayEmbeddings
using DelimitedFiles
using StatsBase
using HypothesisTests

include("../../../src/pecora_uzal_method.jl")
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

using Plots

p1 = histogram(RQA_tde1[1,:],title = "TDE 1 ENTR")
p2 = histogram(RQA_tde1[1,:],title = "GA 1 ENTR")
p3 = histogram(RQA_tde1[1,:],title = "MDOP 1 ENTR")
p4 = histogram(RQA_tde1[1,:],title = "PEC 1 ENTR")
p5 = histogram(RQA_tde2[1,:],title = "TDE 2 ENTR")
p6 = histogram(RQA_tde2[1,:],title = "GA 2 ENTR")
p7 = histogram(RQA_tde2[1,:],title = "MDOP 2 ENTR")
p8 = histogram(RQA_tde2[1,:],title = "PEC 2 ENTR")
plot(p1,p2,p3,p4,p5,p6,p7,p8,layout = (2, 4), legend = false)

p1 = histogram(RQA_tde1[2,:],title = "TDE 1 TRANS")
p2 = histogram(RQA_tde1[2,:],title = "GA 1 TRANS")
p3 = histogram(RQA_tde1[2,:],title = "MDOP 1 TRANS")
p4 = histogram(RQA_tde1[2,:],title = "PEC 1 TRANS")
p5 = histogram(RQA_tde2[2,:],title = "TDE 2 TRANS")
p6 = histogram(RQA_tde2[2,:],title = "GA 2 TRANS")
p7 = histogram(RQA_tde2[2,:],title = "MDOP 2 TRANS")
p8 = histogram(RQA_tde2[2,:],title = "PEC 2 TRANS")
plot(p1,p2,p3,p4,p5,p6,p7,p8,layout = (2, 4), legend = false)

med_TDE = zeros(2,4)
med_GA = zeros(2,4)
med_mdop = zeros(2,4)
med_pec = zeros(2,4)
std_TDE = zeros(2,4)
std_GA = zeros(2,4)
std_mdop = zeros(2,4)
std_pec = zeros(2,4)

for i = 1:4
    # med_TDE[1,i] = median(RQA_tde1[i,:])
    # med_TDE[2,i] = median(RQA_tde2[i,:])
    # med_GA[1,i] = median(RQA_GA1[i,:])
    # med_GA[2,i] = median(RQA_GA2[i,:])
    # med_mdop[1,i] = median(RQA_mdop1[i,:])
    # med_mdop[2,i] = median(RQA_mdop2[i,:])
    # med_pec[1,i] = median(RQA_pec1[i,:])
    # med_pec[2,i] = median(RQA_pec2[i,:])

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

# 1: ENTR
# 2: LAM
# 3: RTE
# 4: TRANS

wtest_tde = UnequalVarianceTTest(RQA_tde1[4,:], RQA_tde2[4,:])
wtest_GA = UnequalVarianceTTest(RQA_GA1[4,:], RQA_GA2[4,:])
wtest_mdop = UnequalVarianceTTest(RQA_mdop1[4,:], RQA_mdop2[4,:])
wtest_pec = UnequalVarianceTTest(RQA_pec1[4,:], RQA_pec2[4,:])
