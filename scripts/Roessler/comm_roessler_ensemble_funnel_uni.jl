using DifferentialEquations
using DynamicalSystems
using DelayEmbeddings
using DelimitedFiles
using Random

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

# Here we look at a chaotic system: The Roessler attractor in the funnel regime

## Integrate system and determine theiler window
Random.seed!(234)
# set time interval for integration
N = 10000 # number of samples
transients = 2000
dt_tra = 0.03
tspan = (dt_tra*N) + (dt_tra*transients)

function roessler!(du,u,p,t)
    x,y,z = u
    a,b,c = p
    du[1] = -(y+z) # dx
    du[2] = x + (a*y) # dy
    du[3] = b+z*(x-c) # dz
end

function prob_func_roe(prob,i,repeat)
  # remake(prob, u0 = [2*rand(1); 2*rand(1); 2*rand(1)])
  # prob
  @. prob.u0 = [2*rand(1); 2*rand(1); 2*rand(1)]
  prob
end

# parameter values
a = .2925
b = .1
c = 8.5
p = [a, b, c]

u0 = [2*rand(1); 2*rand(1); 2*rand(1)]
roe = ODEProblem(roessler!,u0,tspan,p)

roe_ensemble = EnsembleProblem(roe; prob_func = prob_func_roe)
sim = solve(roe_ensemble; saveat = dt_tra, trajectories = 1000)


## Analysis of ensemble

# noise levels
σs = [0, .1]

#preallocation
L_ref = zeros(length(σs),1000)
L_tde = zeros(length(σs),1000)
L_mdop = zeros(length(σs),1000)
L_GA = zeros(length(σs),1000)
L_pec = zeros(length(σs),1000)

dim_tde = zeros(length(σs),1000)
dim_mdop = zeros(length(σs),1000)
dim_GA = zeros(length(σs),1000)
dim_pec = zeros(length(σs),1000)

mfnn_tde = zeros(length(σs),1000)
mfnn_mdop = zeros(length(σs),1000)
mfnn_GA = zeros(length(σs),1000)
mfnn_pec = zeros(length(σs),1000)

jrrf_tde = zeros(length(σs),1000)
jrrf_mdop = zeros(length(σs),1000)
jrrf_GA = zeros(length(σs),1000)
jrrf_pec = zeros(length(σs),1000)

DET_ref = zeros(length(σs),1000)
DET_tde = zeros(length(σs),1000)
DET_mdop = zeros(length(σs),1000)
DET_GA = zeros(length(σs),1000)
DET_pec = zeros(length(σs),1000)

ENTR_ref = zeros(length(σs),1000)
ENTR_tde = zeros(length(σs),1000)
ENTR_mdop = zeros(length(σs),1000)
ENTR_GA = zeros(length(σs),1000)
ENTR_pec = zeros(length(σs),1000)

RTE_ref = zeros(length(σs),1000)
RTE_tde = zeros(length(σs),1000)
RTE_mdop = zeros(length(σs),1000)
RTE_GA = zeros(length(σs),1000)
RTE_pec = zeros(length(σs),1000)

LAM_ref = zeros(length(σs),1000)
LAM_tde = zeros(length(σs),1000)
LAM_mdop = zeros(length(σs),1000)
LAM_GA = zeros(length(σs),1000)
LAM_pec = zeros(length(σs),1000)

TRANS_ref = zeros(length(σs),1000)
TRANS_tde = zeros(length(σs),1000)
TRANS_mdop = zeros(length(σs),1000)
TRANS_GA = zeros(length(σs),1000)
TRANS_pec = zeros(length(σs),1000)

@time for (cnt,σ) in enumerate(σs)
    for i = 1:1000
        display("σ=$σ, i=$i")
        tr = Array(sim[i])'[transients+1:end,:].+σ*randn(10001,3)
        tr = regularize(Dataset(tr))

        if i == 1
            w2 = estimate_delay(tr[:,2], "mi_min")
            global w = w2
            global taus = 0:100

            tw2 = mdop_maximum_delay(tr[:,2])
            lm2 = DelayEmbeddings.findlocalminima(tw2[2])
            if lm2[1] ≤ 2
                lmm2 = try lm2[2]
                catch
                    lm2[1]
                end
            else
                lmm2 = lm2[1]
            end
            global tau_max = lmm2
            global taus_mdop = 0:tau_max
        end

        ## Perform reconstructions and compute according L-value

        #Standard TDE
        Y_tde, τ_tde, _ = optimal_traditional_de(tr[:,2]; τs = taus)
        dim_tde[cnt,i] = size(Y_tde,2)
        τ_tdes = [(i-1)*τ_tde for i = 1:(size(Y_tde,2))]
        L_tde[cnt,i] = compute_delta_L(tr[:,2], τ_tdes, taus[end]; w = w)

        #MDOP
        Y_mdop, τ_vals_mdop, ts_vals_mdop, FNNs_mdop , βs_mdop = mdop_embedding(tr[:,2];
                                                            τs = taus_mdop , w = w)
        dim_mdop[cnt,i] = size(Y_mdop,2)
        L_mdop[cnt,i] = compute_delta_L(tr[:,2], τ_vals_mdop, taus_mdop[end]; w = w)

        #Garcia & Almeida
        Y_GA, τ_vals_GA, ts_vals_GA, FNNs_GA , ns_GA = garcia_almeida_embedding(tr[:,2];
                                                            τs = taus , w = w, T = w)
        dim_GA[cnt,i] = size(Y_GA,2)
        L_GA[cnt,i] = compute_delta_L(tr[:,2], τ_vals_GA, taus[end]; w = w)

        #Pecuzal
        Y_pec, τ_vals_pec, ts_vals_pec, Ls_pec , εs_pec = pecuzal_embedding_update(tr[:,2];
                                                                    τs = taus , w = w)
        dim_pec[cnt,i] = size(Y_pec,2)
        L_pec[cnt,i] = sum(Ls_pec)

        ## compute other evaluation statistics

        mfnn_tde[cnt,i], mfnn_mdop[cnt,i], mfnn_GA[cnt,i], mfnn_pec[cnt,i],
        jrrf_tde[cnt,i], jrrf_mdop[cnt,i], jrrf_GA[cnt,i], jrrf_pec[cnt,i],
        RQA_ref, RQA1, RQA2, RQA3, RQA4, _, _, _, _, _ =
                perform_recurrence_analysis(tr, Dataset(Y_tde), Dataset(Y_mdop),
                            Dataset(Y_GA), Dataset(Y_pec); ε = 0.08, w = w, kNN = 10)


        DET_ref[cnt,i] = RQA_ref.DET
        RTE_ref[cnt,i] = RQA_ref.RTE
        ENTR_ref[cnt,i] = RQA_ref.ENTR
        LAM_ref[cnt,i] = RQA_ref.LAM
        TRANS_ref[cnt,i] = RQA_ref.TRANS

        DET_tde[cnt,i] = RQA1.DET
        RTE_tde[cnt,i] = RQA1.RTE
        ENTR_tde[cnt,i] = RQA1.ENTR
        LAM_tde[cnt,i] = RQA1.LAM
        TRANS_tde[cnt,i] = RQA1.TRANS

        DET_mdop[cnt,i] = RQA2.DET
        RTE_mdop[cnt,i] = RQA2.RTE
        ENTR_mdop[cnt,i] = RQA2.ENTR
        LAM_mdop[cnt,i] = RQA2.LAM
        TRANS_mdop[cnt,i] = RQA2.TRANS

        DET_GA[cnt,i] = RQA3.DET
        RTE_GA[cnt,i] = RQA3.RTE
        ENTR_GA[cnt,i] = RQA3.ENTR
        LAM_GA[cnt,i] = RQA3.LAM
        TRANS_GA[cnt,i] = RQA3.TRANS

        DET_pec[cnt,i] = RQA4.DET
        RTE_pec[cnt,i] = RQA4.RTE
        ENTR_pec[cnt,i] = RQA4.ENTR
        LAM_pec[cnt,i] = RQA4.LAM
        TRANS_pec[cnt,i] = RQA4.TRANS
    end
end

# save results
writedlm("./ensemble results funnel uni/L_tde.csv", L_tde)
writedlm("./ensemble results funnel uni/L_mdop.csv", L_mdop)
writedlm("./ensemble results funnel uni/L_GA.csv", L_GA)
writedlm("./ensemble results funnel uni/L_pec.csv", L_pec)

writedlm("./ensemble results funnel uni/dim_tde.csv", dim_tde)
writedlm("./ensemble results funnel uni/dim_mdop.csv", dim_mdop)
writedlm("./ensemble results funnel uni/dim_GA.csv", dim_GA)
writedlm("./ensemble results funnel uni/dim_pec.csv", dim_pec)

writedlm("./ensemble results funnel uni/mfnn_tde.csv", mfnn_tde)
writedlm("./ensemble results funnel uni/mfnn_mdop.csv", mfnn_mdop)
writedlm("./ensemble results funnel uni/mfnn_GA.csv", mfnn_GA)
writedlm("./ensemble results funnel uni/mfnn_pec.csv", mfnn_pec)

writedlm("./ensemble results funnel uni/jrrf_tde.csv", jrrf_tde)
writedlm("./ensemble results funnel uni/jrrf_mdop.csv", jrrf_mdop)
writedlm("./ensemble results funnel uni/jrrf_GA.csv", jrrf_GA)
writedlm("./ensemble results funnel uni/jrrf_pec.csv", jrrf_pec)

writedlm("./ensemble results funnel uni/DET_ref.csv", DET_ref)
writedlm("./ensemble results funnel uni/DET_tde.csv", DET_tde)
writedlm("./ensemble results funnel uni/DET_mdop.csv", DET_mdop)
writedlm("./ensemble results funnel uni/DET_GA.csv", DET_GA)
writedlm("./ensemble results funnel uni/DET_pec.csv", DET_pec)

writedlm("./ensemble results funnel uni/ENTR_ref.csv", ENTR_ref)
writedlm("./ensemble results funnel uni/ENTR_tde.csv", ENTR_tde)
writedlm("./ensemble results funnel uni/ENTR_mdop.csv", ENTR_mdop)
writedlm("./ensemble results funnel uni/ENTR_GA.csv", ENTR_GA)
writedlm("./ensemble results funnel uni/ENTR_pec.csv", ENTR_pec)

writedlm("./ensemble results funnel uni/RTE_ref.csv", RTE_ref)
writedlm("./ensemble results funnel uni/RTE_tde.csv", RTE_tde)
writedlm("./ensemble results funnel uni/RTE_mdop.csv", RTE_mdop)
writedlm("./ensemble results funnel uni/RTE_GA.csv", RTE_GA)
writedlm("./ensemble results funnel uni/RTE_pec.csv", RTE_pec)

writedlm("./ensemble results funnel uni/LAM_ref.csv", LAM_ref)
writedlm("./ensemble results funnel uni/LAM_tde.csv", LAM_tde)
writedlm("./ensemble results funnel uni/LAM_mdop.csv", LAM_mdop)
writedlm("./ensemble results funnel uni/LAM_GA.csv", LAM_GA)
writedlm("./ensemble results funnel uni/LAM_pec.csv", LAM_pec)

writedlm("./ensemble results funnel uni/TRANS_ref.csv", TRANS_ref)
writedlm("./ensemble results funnel uni/TRANS_tde.csv", TRANS_tde)
writedlm("./ensemble results funnel uni/TRANS_mdop.csv", TRANS_mdop)
writedlm("./ensemble results funnel uni/TRANS_GA.csv", TRANS_GA)
writedlm("./ensemble results funnel uni/TRANS_pec.csv", TRANS_pec)
