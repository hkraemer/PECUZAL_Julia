using DrWatson
@quickactivate "new-embedding-methods"

using DifferentialEquations
using DynamicalSystems
using DelayEmbeddings
using DelimitedFiles

include("../../src/pecora_uzal_method.jl")
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

# Here we look at a chaotic system: The Mackey Glass system

## 

# set time interval for integration
N = 10000 # number of samples
transients = 2000
dt_tra = 0.5
tspan = (dt_tra*N) + (dt_tra*transients)


function mackey_glass(du,u,h,p,t)
  beta,n,gamma,tau = p
  hist = h(p, t-tau)[1]
  du[1] = (beta*hist)/(1+hist^n) - gamma * u[1]
end

function prob_func_mackey(prob,i,repeat)
  @. prob.u0 = [1.5*rand(1)]
  prob
end

# set parameters
h(p,t) = 0.
tau_d = 44.
n = 10.
β = 0.2
γ = 0.1

p = (β,n,γ,tau_d)

u0 = [1.0]


alg = MethodOfSteps(Tsit5())

## Analysis of ensemble

# method for standard τ-estimation
method = "mi_min"
#method = "ac_zero"

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

@time for (cnt,σ) in enumerate(σs)
    for i = 1:1000
        display("σ=$σ, i=$i")
        mackey = DDEProblem(mackey_glass,u0,h,tspan,p; constant_lags=tau_d)
        sol = solve(mackey,alg; adaptive=false, dt=dt_tra)

        tr = Array(sol)'[transients+1:end,:].+σ*randn(10001,3)
        tr = regularize(Dataset(tr))

        if i == 1
            w1 = estimate_delay(tr[:,1], method)
            global w = w1
            global taus = 0:100

            tw1 = mdop_maximum_delay(tr[:,1])
            lm1 = DelayEmbeddings.findlocalminima(tw1[2])
            if lm1[1] ≤ 2
                lmm1 = try lm1[2]
                catch
                    lm1[1]
                end
            else
                lmm1 = lm1[1]
            end

            global tau_max = lmm1
            global taus_mdop = 0:tau_max
        end

        ## Perform reconstructions and compute according L-value

        # L-value of reference
        L_ref[cnt,i] = uzal_cost(tr, Tw = (4*w), w = w, samplesize=1)

        #Standard TDE
        Y_tde, τ_tde = standard_embedding_cao(tr[:,1], method = method)
        dim_tde[cnt,i] = size(Y_tde,2)
        L_tde[cnt,i] = uzal_cost(Y_tde, Tw = (4*w), w = w, samplesize=1)

        #MDOP
        Y_mdop, τ_vals_mdop, ts_vals_mdop, FNNs_mdop , βs_mdop = mdop_embedding(tr[:,1];
                                                            τs = taus_mdop , w = w)
        dim_mdop[cnt,i] = size(Y_mdop,2)
        L_mdop[cnt,i] = uzal_cost(Y_mdop, Tw = (4*w), w = w, samplesize=1)

        #Garcia & Almeida
        Y_GA, τ_vals_GA, ts_vals_GA, FNNs_GA , ns_GA = garcia_almeida_embedding(tr[:,1];
                                                            τs = taus , w = w, T = w)
        dim_GA[cnt,i] = size(Y_GA,2)
        L_GA[cnt,i] = uzal_cost(Y_GA, Tw = (4*w), w = w, samplesize=1)

        #Pecuzal
        Y_pec, τ_vals_pec, ts_vals_pec, Ls_pec , εs_pec = pecuzal_embedding(tr[:,1];
                                                                    τs = taus , w = w)
        dim_pec[cnt,i] = size(Y_pec,2)
        L_pec[cnt,i] = minimum(Ls_pec)

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

        DET_tde[cnt,i] = RQA1.DET
        RTE_tde[cnt,i] = RQA1.RTE
        ENTR_tde[cnt,i] = RQA1.ENTR
        LAM_tde[cnt,i] = RQA1.LAM

        DET_mdop[cnt,i] = RQA2.DET
        RTE_mdop[cnt,i] = RQA2.RTE
        ENTR_mdop[cnt,i] = RQA2.ENTR
        LAM_mdop[cnt,i] = RQA2.LAM

        DET_GA[cnt,i] = RQA3.DET
        RTE_GA[cnt,i] = RQA3.RTE
        ENTR_GA[cnt,i] = RQA3.ENTR
        LAM_GA[cnt,i] = RQA3.LAM

        DET_pec[cnt,i] = RQA4.DET
        RTE_pec[cnt,i] = RQA4.RTE
        ENTR_pec[cnt,i] = RQA4.ENTR
        LAM_pec[cnt,i] = RQA4.LAM
    end
end

# save results
writedlm("./scripts/Mackey Glass/ensemble results/L_ref.csv", L_ref)
writedlm("./scripts/Mackey Glass/ensemble results/L_tde.csv", L_tde)
writedlm("./scripts/Mackey Glass/ensemble results/L_mdop.csv", L_mdop)
writedlm("./scripts/Mackey Glass/ensemble results/L_GA.csv", L_GA)
writedlm("./scripts/Mackey Glass/ensemble results/L_pec.csv", L_pec)

writedlm("./scripts/Mackey Glass/ensemble results/dim_tde.csv", dim_tde)
writedlm("./scripts/Mackey Glass/ensemble results/dim_mdop.csv", dim_mdop)
writedlm("./scripts/Mackey Glass/ensemble results/dim_GA.csv", dim_GA)
writedlm("./scripts/Mackey Glass/ensemble results/dim_pec.csv", dim_pec)

writedlm("./scripts/Mackey Glass/ensemble results/mfnn_tde.csv", mfnn_tde)
writedlm("./scripts/Mackey Glass/ensemble results/mfnn_mdop.csv", mfnn_mdop)
writedlm("./scripts/Mackey Glass/ensemble results/mfnn_GA.csv", mfnn_GA)
writedlm("./scripts/Mackey Glass/ensemble results/mfnn_pec.csv", mfnn_pec)

writedlm("./scripts/Mackey Glass/ensemble results/jrrf_tde.csv", jrrf_tde)
writedlm("./scripts/Mackey Glass/ensemble results/jrrf_mdop.csv", jrrf_mdop)
writedlm("./scripts/Mackey Glass/ensemble results/jrrf_GA.csv", jrrf_GA)
writedlm("./scripts/Mackey Glass/ensemble results/jrrf_pec.csv", jrrf_pec)

writedlm("./scripts/Mackey Glass/ensemble results/DET_ref.csv", DET_ref)
writedlm("./scripts/Mackey Glass/ensemble results/DET_tde.csv", DET_tde)
writedlm("./scripts/Mackey Glass/ensemble results/DET_mdop.csv", DET_mdop)
writedlm("./scripts/Mackey Glass/ensemble results/DET_GA.csv", DET_GA)
writedlm("./scripts/Mackey Glass/ensemble results/DET_pec.csv", DET_pec)

writedlm("./scripts/Mackey Glass/ensemble results/ENTR_ref.csv", ENTR_ref)
writedlm("./scripts/Mackey Glass/ensemble results/ENTR_tde.csv", ENTR_tde)
writedlm("./scripts/Mackey Glass/ensemble results/ENTR_mdop.csv", ENTR_mdop)
writedlm("./scripts/Mackey Glass/ensemble results/ENTR_GA.csv", ENTR_GA)
writedlm("./scripts/Mackey Glass/ensemble results/ENTR_pec.csv", ENTR_pec)

writedlm("./scripts/Mackey Glass/ensemble results/RTE_ref.csv", RTE_ref)
writedlm("./scripts/Mackey Glass/ensemble results/RTE_tde.csv", RTE_tde)
writedlm("./scripts/Mackey Glass/ensemble results/RTE_mdop.csv", RTE_mdop)
writedlm("./scripts/Mackey Glass/ensemble results/RTE_GA.csv", RTE_GA)
writedlm("./scripts/Mackey Glass/ensemble results/RTE_pec.csv", RTE_pec)

writedlm("./scripts/Mackey Glass/ensemble results/LAM_ref.csv", LAM_ref)
writedlm("./scripts/Mackey Glass/ensemble results/LAM_tde.csv", LAM_tde)
writedlm("./scripts/Mackey Glass/ensemble results/LAM_mdop.csv", LAM_mdop)
writedlm("./scripts/Mackey Glass/ensemble results/LAM_GA.csv", LAM_GA)
writedlm("./scripts/Mackey Glass/ensemble results/LAM_pec.csv", LAM_pec)
