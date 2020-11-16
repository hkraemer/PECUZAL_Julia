using DrWatson
@quickactivate "new-embedding-methods"

using DelayEmbeddings
using DelimitedFiles
using StatsBase
using PyPlot
pygui(true)

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

data1 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/s3ch4.txt")
data2 = readdlm("./scripts/Experimental data/Standard electrochemical chaos/s1ch3.txt")

# downsampling
x1 = data1[1:2:end]
x2 = data2[1:2:end]

x1 = regularize(x1 .+ .00001*randn(length(x1)))
x2 = regularize(x2 .+ .00001*randn(length(x2)))

w1 = estimate_delay(x1, "mi_min")
w2 = estimate_delay(x2, "mi_min")

taus = 0:150

## Perform embeddings

#Standard TDE
@time Y_tde1, _ = standard_embedding_cao(x1)
@time Y_tde2, _ = standard_embedding_cao(x2)
L_tde1 = uzal_cost(Y_tde1, Tw = (4*w1), w = w1, samplesize=1)
L_tde2 = uzal_cost(Y_tde2, Tw = (4*w2), w = w2, samplesize=1)


#MDOP embedding
println("Computation time MDOP:")
@time begin
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
    L_mdop1 = uzal_cost(Y_mdop1, Tw = (4*w1), w = w1, samplesize=1)
    Y_mdop2, τ_vals_mdop2, ts_vals_mdop2, FNNs_mdop2 , βs_mdop2 = mdop_embedding(x2;
                                                        τs = taus_mdop2 , w = w2)
    L_mdop2 = uzal_cost(Y_mdop2, Tw = (4*w2), w = w2, samplesize=1)
end

#Garcia & Almeida
println("Computation time Garcia & Almeida method:")
@time Y_GA1, τ_vals_GA1, ts_vals_GA1, FNNs_GA1 , ns_GA1 = garcia_almeida_embedding(x1;
                                                    τs = taus , w = w1, T = w1)
L_GA1 = uzal_cost(Y_GA1, Tw = (4*w1), w = w1, samplesize=1)
@time Y_GA2, τ_vals_GA2, ts_vals_GA2, FNNs_GA2 , ns_GA2 = garcia_almeida_embedding(x2;
                                                    τs = taus , w = w2, T = w2)
L_GA2 = uzal_cost(Y_GA2, Tw = (4*w2), w = w2, samplesize=1)

#Pecuzal
println("Computation time pecuzal method:")
@time Y_pec1, τ_vals_pec1, ts_vals_pec1, Ls_pec1 , εs_pec1 = pecuzal_embedding(x1;
                                                            τs = taus , w = w1)
L_pec1 = minimum(Ls_pec1)
@time Y_pec2, τ_vals_pec2, ts_vals_pec2, Ls_pec2 , εs_pec2 = pecuzal_embedding(x2;
                                                            τs = taus , w = w2)
L_pec2 = minimum(Ls_pec2)

## Plot attractor
figure()
plot3D(Y_pec1[:,1],Y_pec1[:,2],Y_pec1[:,3], linewidth=.1)
title("PECUZAL reconstruction of chaotic chemical oscillator (setup I)")
xlabel("x(t+$(τ_vals_pec1[1]))")
ylabel("x(t+$(τ_vals_pec1[2]))")
zlabel("x(t+$(τ_vals_pec1[3]))")
grid()

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

# set recurrence threshold for further analysis
ε = 0.08

# set length of time series sample
N = 2000

## RQA
samplesize = 1000
idxs1 = sample(vec(1:NN1-N-1), samplesize, replace = false)
idxs2 = sample(vec(1:NN2-N-1), samplesize, replace = false)

lmin = 2

RQA_tde1 = zeros(4,samplesize)
RQA_tde2 = zeros(4,samplesize)
RQA_GA1 = zeros(4,samplesize)
RQA_GA2 = zeros(4,samplesize)
RQA_mdop1 = zeros(4,samplesize)
RQA_mdop2 = zeros(4,samplesize)
RQA_pec1 = zeros(4,samplesize)
RQA_pec2 = zeros(4,samplesize)

for i = 1:samplesize
    display("run: $i")

    idx1 = idxs1[i]
    idx2 = idxs2[i]

    R11 = RecurrenceMatrix(Y_tde1[idx1:idx1+N,:], ε; fixedrate = true)
    R21 = RecurrenceMatrix(Y_GA1[idx1:idx1+N,:], ε; fixedrate = true)
    R31 = RecurrenceMatrix(Y_mdop1[idx1:idx1+N,:], ε; fixedrate = true)
    R41 = RecurrenceMatrix(Y_pec1[idx1:idx1+N,:], ε; fixedrate = true)

    R12 = RecurrenceMatrix(Y_tde2[idx2:idx2+N,:], ε; fixedrate = true)
    R22 = RecurrenceMatrix(Y_GA2[idx2:idx2+N,:], ε; fixedrate = true)
    R32 = RecurrenceMatrix(Y_mdop2[idx2:idx2+N,:], ε; fixedrate = true)
    R42 = RecurrenceMatrix(Y_pec2[idx2:idx2+N,:], ε; fixedrate = true)

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
    RQA_tde1[4,i] = transitivity(R11)
    RQA_tde2[1,i] = RQA11.ENTR
    RQA_tde2[2,i] = RQA11.LAM
    RQA_tde2[3,i] = RQA11.RTE
    RQA_tde2[4,i] = transitivity(R12)

    RQA_GA1[1,i] = RQA21.ENTR
    RQA_GA1[2,i] = RQA21.LAM
    RQA_GA1[3,i] = RQA21.RTE
    RQA_GA1[4,i] = transitivity(R21)
    RQA_GA2[1,i] = RQA22.ENTR
    RQA_GA2[2,i] = RQA22.LAM
    RQA_GA2[3,i] = RQA22.RTE
    RQA_GA2[4,i] = transitivity(R22)

    RQA_mdop1[1,i] = RQA31.ENTR
    RQA_mdop1[2,i] = RQA31.LAM
    RQA_mdop1[3,i] = RQA31.RTE
    RQA_mdop1[4,i] = transitivity(R31)
    RQA_mdop2[1,i] = RQA32.ENTR
    RQA_mdop2[2,i] = RQA32.LAM
    RQA_mdop2[3,i] = RQA32.RTE
    RQA_mdop2[4,i] = transitivity(R32)

    RQA_pec1[1,i] = RQA41.ENTR
    RQA_pec1[2,i] = RQA41.LAM
    RQA_pec1[3,i] = RQA41.RTE
    RQA_pec1[4,i] = transitivity(R41)
    RQA_pec2[1,i] = RQA42.ENTR
    RQA_pec2[2,i] = RQA42.LAM
    RQA_pec2[3,i] = RQA42.RTE
    RQA_pec2[4,i] = transitivity(R42)
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


##
