using DrWatson
@quickactivate "PECUZAL_Julia"

using DelayEmbeddings
using DelimitedFiles
using StatsBase
using Random

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

## Load the data
Random.seed!(1234)

data1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/s3ch4.txt")
data2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/s1ch3.txt")

# length of the whole time series
M1 = length(data1)
M2 = length(data2)
M = minimum(hcat(M1, M2))

# downsampling
data1 = data1[1:4:M]
data2 = data2[1:4:M]
M = length(data1)

data1 = regularize(data1 .+ .000001*randn(length(data1)))
data2 = regularize(data2 .+ .000001*randn(length(data2)))

# the possible lags
taus = 0:100

# set recurrence threshold for further analysis (and min line length)
ε = 0.08
lmin = 2

# set length of time series sample
N = 1000

# number of time series samples
samplesize = 1000

# draw samples
idxs1 = sample(vec(1:M-N-1), samplesize, replace = false)
idxs2 = sample(vec(1:M-N-1), samplesize, replace = false)

#preallocation
RQA_tde1 = zeros(4,samplesize)
RQA_tde2 = zeros(4,samplesize)
RQA_GA1 = zeros(4,samplesize)
RQA_GA2 = zeros(4,samplesize)
RQA_mdop1 = zeros(4,samplesize)
RQA_mdop2 = zeros(4,samplesize)
RQA_pec1 = zeros(4,samplesize)
RQA_pec2 = zeros(4,samplesize)

dim_tde1 = zeros(Int, samplesize)
dim_tde2 = zeros(Int, samplesize)
L_tde1 = zeros(samplesize)
L_tde2 = zeros(samplesize)

dim_GA1 = zeros(Int, samplesize)
dim_GA2 = zeros(Int, samplesize)
L_GA1 = zeros(samplesize)
L_GA2 = zeros(samplesize)

dim_mdop1 = zeros(Int, samplesize)
dim_mdop2 = zeros(Int, samplesize)
L_mdop1 = zeros(samplesize)
L_mdop2 = zeros(samplesize)

dim_pec1 = zeros(Int, samplesize)
dim_pec2 = zeros(Int, samplesize)
L_pec1 = zeros(samplesize)
L_pec2 = zeros(samplesize)

for i = 1:samplesize
    display("run: $i")
    # draw the actual sample
    idx1 = idxs1[i]
    idx2 = idxs2[i]
    x1 = data1[idx1:idx1+N]
    x2 = data2[idx2:idx2+N]

    w1 = estimate_delay(x1, "mi_min")
    w2 = estimate_delay(x2, "mi_min")

    ## Perform embeddings
    #Standard TDE
    Y_tde1, τ_tde1, _ = optimal_traditional_de(x1; τs = taus)
    Y_tde2, τ_tde2, _ = optimal_traditional_de(x2; τs = taus)

    dim_tde1[i] = size(Y_tde1,2)
    dim_tde2[i] = size(Y_tde2,2)
    τ_tdes1 = [(i-1)*τ_tde1 for i = 1:(size(Y_tde1,2))]
    τ_tdes2 = [(i-1)*τ_tde2 for i = 1:(size(Y_tde2,2))]
    L_tde1[i] = compute_delta_L(x1, τ_tdes1, taus[end]; w = w1)
    L_tde2[i] = compute_delta_L(x2, τ_tdes2, taus[end]; w = w2)

    #MDOP embedding
    tw1 = mdop_maximum_delay(x1)
    tw2 = mdop_maximum_delay(x2)
    lm1 = DelayEmbeddings.findlocalminima(tw1[2])
    if lm1[1] ≤ 3
        lmm1 = try lm1[2]
        catch
            lm1[1]
        end
    else
        lmm1 = lm1[1]
    end
    lm2 = DelayEmbeddings.findlocalminima(tw2[2])
    if lm2[1] ≤ 3
        lmm2 = try lm2[2]
        catch
            lm2[1]
        end
    else
        lmm2 = lm2[1]
    end
    tau_max1 = lmm1
    tau_max2 = lmm2
    taus_mdop1 = 0:tau_max1
    taus_mdop2 = 0:tau_max2

    Y_mdop1, τ_vals_mdop1, ts_vals_mdop1, FNNs_mdop1 , βs_mdop1 = mdop_embedding(x1;
                                                        τs = taus_mdop1 , w = w1)
    dim_mdop1[i] = size(Y_mdop1,2)
    L_mdop1[i] = compute_delta_L(x1, τ_vals_mdop1, taus_mdop1[end]; w = w1)

    Y_mdop2, τ_vals_mdop2, ts_vals_mdop2, FNNs_mdop2 , βs_mdop2 = mdop_embedding(x2;
                                                        τs = taus_mdop2 , w = w2)
    dim_mdop2[i] = size(Y_mdop2,2)
    L_mdop2[i] = compute_delta_L(x2, τ_vals_mdop2, taus_mdop2[end]; w = w2)

    #Garcia & Almeida

    Y_GA1, τ_vals_GA1, ts_vals_GA1, FNNs_GA1 , ns_GA1 = garcia_almeida_embedding(x1;
                                                        τs = taus , w = w1, T = w1)
    dim_GA1[i] = size(Y_GA1,2)
    L_GA1[i] = compute_delta_L(x1, τ_vals_GA1, taus[end]; w = w1)
    Y_GA2, τ_vals_GA2, ts_vals_GA2, FNNs_GA2 , ns_GA2 = garcia_almeida_embedding(x2;
                                                        τs = taus , w = w2, T = w2)
    dim_GA2[i] = size(Y_GA2,2)
    L_GA2[i] = compute_delta_L(x2, τ_vals_GA2, taus[end]; w = w2)

    #Pecuzal
    Y_pec1, τ_vals_pec1, ts_vals_pec1, Ls_pec1 , εs_pec1 = pecuzal_embedding_update(x1;
                                                                τs = taus , w = w1)
    dim_pec1[i] = size(Y_pec1,2)
    L_pec1[i] = sum(Ls_pec1)
    Y_pec2, τ_vals_pec2, ts_vals_pec2, Ls_pec2 , εs_pec2 = pecuzal_embedding_update(x2;
                                                                τs = taus , w = w2)
    dim_pec2[i] = size(Y_pec2,2)
    L_pec2[i] = sum(Ls_pec2)

    ## Make all reconstructions equally long
    N1 = length(Y_tde1)
    N2 = length(Y_GA1)
    N3 = length(Y_mdop1)
    N4 = length(Y_pec1)
    NN1 = minimum(hcat(N1, N2, N3, N4))

    N1 = length(Y_tde2)
    N2 = length(Y_GA2)
    N3 = length(Y_mdop2)
    N4 = length(Y_pec2)
    NN2 = minimum(hcat(N1, N2, N3, N4))

    ## RQA
    R11 = RecurrenceMatrix(Y_tde1[1:NN1,:], ε; fixedrate = true)
    R21 = RecurrenceMatrix(Y_GA1[1:NN1,:], ε; fixedrate = true)
    R31 = RecurrenceMatrix(Y_mdop1[1:NN1,:], ε; fixedrate = true)
    R41 = RecurrenceMatrix(Y_pec1[1:NN1,:], ε; fixedrate = true)

    R12 = RecurrenceMatrix(Y_tde2[1:NN2,:], ε; fixedrate = true)
    R22 = RecurrenceMatrix(Y_GA2[1:NN2,:], ε; fixedrate = true)
    R32 = RecurrenceMatrix(Y_mdop2[1:NN2,:], ε; fixedrate = true)
    R42 = RecurrenceMatrix(Y_pec2[1:NN2,:], ε; fixedrate = true)

    RQA11 = rqa(R11; theiler = w1, lmin = lmin)
    RQA21 = rqa(R21; theiler = w1, lmin = lmin)
    RQA31 = rqa(R31; theiler = w1, lmin = lmin)
    RQA41 = rqa(R41; theiler = w1, lmin = lmin)

    RQA12 = rqa(R12; theiler = w2, lmin = lmin)
    RQA22 = rqa(R22; theiler = w2, lmin = lmin)
    RQA32 = rqa(R32; theiler = w2, lmin = lmin)
    RQA42 = rqa(R42; theiler = w2, lmin = lmin)

    RQA_tde1[1,i] = RQA11.ENTR
    RQA_tde1[2,i] = RQA11.LAM
    RQA_tde1[3,i] = RQA11.RTE
    RQA_tde1[4,i] = RQA11.TRANS
    RQA_tde2[1,i] = RQA12.ENTR
    RQA_tde2[2,i] = RQA12.LAM
    RQA_tde2[3,i] = RQA12.RTE
    RQA_tde2[4,i] = RQA12.TRANS

    RQA_GA1[1,i] = RQA21.ENTR
    RQA_GA1[2,i] = RQA21.LAM
    RQA_GA1[3,i] = RQA21.RTE
    RQA_GA1[4,i] = RQA21.TRANS
    RQA_GA2[1,i] = RQA22.ENTR
    RQA_GA2[2,i] = RQA22.LAM
    RQA_GA2[3,i] = RQA22.RTE
    RQA_GA2[4,i] = RQA22.TRANS

    RQA_mdop1[1,i] = RQA31.ENTR
    RQA_mdop1[2,i] = RQA31.LAM
    RQA_mdop1[3,i] = RQA31.RTE
    RQA_mdop1[4,i] = RQA31.TRANS
    RQA_mdop2[1,i] = RQA32.ENTR
    RQA_mdop2[2,i] = RQA32.LAM
    RQA_mdop2[3,i] = RQA32.RTE
    RQA_mdop2[4,i] = RQA32.TRANS

    RQA_pec1[1,i] = RQA41.ENTR
    RQA_pec1[2,i] = RQA41.LAM
    RQA_pec1[3,i] = RQA41.RTE
    RQA_pec1[4,i] = RQA41.TRANS
    RQA_pec2[1,i] = RQA42.ENTR
    RQA_pec2[2,i] = RQA42.LAM
    RQA_pec2[3,i] = RQA42.RTE
    RQA_pec2[4,i] = RQA42.TRANS
end

# save results
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_tde1.csv", RQA_tde1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_tde2.csv", RQA_tde2)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_GA1.csv", RQA_GA1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_GA2.csv", RQA_GA2)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_mdop1.csv", RQA_mdop1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_mdop2.csv", RQA_mdop2)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_pec1.csv", RQA_pec1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/RQA_pec2.csv", RQA_pec2)

writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_tde1.csv", dim_tde1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_tde2.csv", dim_tde2)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_tde1.csv", L_tde1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_tde2.csv", L_tde2)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_GA1.csv", dim_GA1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_GA2.csv", dim_GA2)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_GA1.csv", L_GA1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_GA2.csv", L_GA2)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_mdop1.csv", dim_mdop1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_mdop2.csv", dim_mdop2)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_mdop1.csv", L_mdop1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_mdop2.csv", L_mdop2)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_pec1.csv", dim_pec1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/dim_pec2.csv", dim_pec2)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_pec1.csv", L_pec1)
writedlm("./scripts/Experimental data/Standard electrochemical chaos/RQA results/L_pec2.csv", L_pec2)

##
# println("TDE:")
# println(mean(dim_tde1))
# println(mean(dim_tde2))
# println(mean(L_tde1))
# println(mean(L_tde2))
# println("***********")
# println("GA:")
# println(mean(dim_GA1))
# println(mean(dim_GA2))
# println(mean(L_GA1))
# println(mean(L_GA2))
# println("***********")
# println("mdop:")
# println(mean(dim_mdop1))
# println(mean(dim_mdop2))
# println(mean(L_mdop1))
# println(mean(L_mdop2))
# println("***********")
# println("pec:")
# println(mean(dim_pec1))
# println(mean(dim_pec2))
# println(mean(L_pec1))
# println(mean(L_pec2))
# println("***********")
