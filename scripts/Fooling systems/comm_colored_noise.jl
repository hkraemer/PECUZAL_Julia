using DrWatson
@quickactivate "new-embedding-methods"

using DifferentialEquations
using DynamicalSystems
using DelayEmbeddings
using DelimitedFiles

include("../../src/pecuzal_method.jl")
include("../../src/data_analysis_functions.jl")

## We analyze the reconstruction from standard time delay embedding, the
# MDOP embedding algorithm and the Garcia & Almeida method and compare it to
# our propsed method: The pecuzal method.
# For evaluation of the different reconstructions we consider Uzal's L-statistic,
# the fraction of recurrence rates from the JRP of the reconstruction and the
# reference and only the reference. As a third criterion we look at the
# recurrence time entropy, since it is related to the Kolmogorov-entropy as a
# dynamical invariant. We also look at two other RQA-measures and, moreover, the
# generalized mutual false nearest neighbors.

# Here we try to fool the reconstruction algorithms by feeding in stochastic
# time series.

## Generate data

N = 5000 # time series length
tr = zeros(N,3) # multivariate dataset we will use
alpha, p = .9, 1.
u0 = .2
tr[:,1] = ar_process(u0, alpha, p, N)
tr[:,2] = randn(N)
tr[:,3] = rand(N)

tr = regularize(Dataset(tr))

#method = "mi_min"
method = "ac_zero"

w1 = estimate_delay(tr[:,1], method)
w2 = estimate_delay(tr[:,2], method)
w3 = estimate_delay(tr[:,3], method)

w = maximum(hcat(w1,w2,w3))
taus = 0:100

# L-value of reference
L_ref = uzal_cost(tr, Tw = (4*w), w = w, samplesize=1)

#Standard TDE
@time Y_tde, τ = standard_embedding_cao(tr[:,1])
E2s = DelayEmbeddings.stochastic_indicator(tr[:,1], τ, 1:size(Y_tde,2))
L_tde = uzal_cost(Y_tde, Tw = (4*w), w = w, samplesize=1)

#MDOP embedding
println("Computation time MDOP:")
@time begin
    tw1 = mdop_maximum_delay(tr[:,1])
    tw2 = mdop_maximum_delay(tr[:,2])
    tw3 = mdop_maximum_delay(tr[:,3])
    lm1 = DelayEmbeddings.findlocalminima(tw1[2])
    if lm1[1] ≤ 2
        lmm1 = try lm1[2]
        catch
            lm1[1]
        end
    else
        lmm1 = lm1[1]
    end
    lm2 = DelayEmbeddings.findlocalminima(tw2[2])
    if lm2[1] ≤ 2
        lmm2 = try lm2[2]
        catch
            lm2[1]
        end
    else
        lmm2 = lm2[1]
    end
    lm3 = DelayEmbeddings.findlocalminima(tw3[2])
    if lm3[1] ≤ 2
        lmm3 = try lm3[2]
        catch
            lm3[1]
        end
    else
        lmm3 = lm3[1]
    end
    tau_max = maximum(hcat(lmm1, lmm2, lmm3))
    taus_mdop = 0:tau_max

    Y_mdop, τ_vals_mdop, ts_vals_mdop, FNNs_mdop , βs_mdop = mdop_embedding(tr;
                                                        τs = taus_mdop , w = w)
    L_mdop = uzal_cost(Y_mdop, Tw = (4*w), w = w, samplesize=1)
end

#Garcia & Almeida
println("Computation time Garcia & Almeida method:")
@time Y_GA, τ_vals_GA, ts_vals_GA, FNNs_GA , ns_GA = garcia_almeida_embedding(tr;
                                                    τs = taus , w = w, T = w)
L_GA = uzal_cost(Y_GA, Tw = (4*w), w = w, samplesize=1)

#Pecuzal
println("Computation time pecuzal method:")
@time Y_pec, τ_vals_pec, ts_vals_pec, Ls_pec , εs_pec = pecuzal_embedding(tr;
                                                            τs = taus , w = w)
L_pec = minimum(Ls_pec)


## compute other evaluation statistics

@time mfnn1, mfnn2, mfnn3, mfnn4, f1, f2, f3, f4, RQA_ref, RQA1, RQA2,
                                            RQA3, RQA4, R_ref, R1, R2, R3, R4 =
        perform_recurrence_analysis(tr, Dataset(Y_tde), Dataset(Y_mdop),
                    Dataset(Y_GA), Dataset(Y_pec); ε = 0.08, w = w, kNN = 10)


# writedlm("./scripts/Fooling systems/noise results/Y_GA.csv",Y_GA)
# writedlm("./scripts/Fooling systems/noise results/Y_mdop.csv",Y_mdop)
# writedlm("./scripts/Fooling systems/noise results/Y_pec.csv",Y_pec)
# writedlm("./scripts/Fooling systems/noise results/taus_GA.csv", τ_vals_GA)
# writedlm("./scripts/Fooling systems/noise results/taus_mdop.csv", τ_vals_mdop)
# writedlm("./scripts/Fooling systems/noise results/taus_pec.csv", τ_vals_pec)
# writedlm("./scripts/Fooling systems/noise results/ts_GA.csv", ts_vals_GA)
# writedlm("./scripts/Fooling systems/noise results/ts_mdop.csv", ts_vals_mdop)
# writedlm("./scripts/Fooling systems/noise results/ts_pec.csv", ts_vals_pec)
