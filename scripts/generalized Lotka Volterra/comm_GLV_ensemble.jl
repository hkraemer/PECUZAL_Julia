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

# Here we look at a complex system: The generalized Lotka-Volterra-system with 9
# specicies.

##

# set time interval for integration
N = 5000 # number of samples
transients = 300
dt_tra = 0.15
tspan = (dt_tra*N) + (dt_tra*transients)


function glv!(du,u,p,t)
    ρ,γ,r,e,f,c,d = p
    m0 = [0 c e; e 0 c; c e 0]
    md = [d r r; r d r; r r d]
    mf = [f r r; r f r; r r f]
    A = [m0 md mf; mf m0 md; md mf m0]
    Amix = zeros(size(A))
    for i = 1:size(A,1)
        for j = 1:size(A,1)
            Amix[i,j] = A[i,j]*u[i]*u[j]
        end
        du[i] = ρ*u[i]-γ*u[i]^2-sum(Amix[i,:])
    end
end

function prob_func_glv(prob,i,repeat)
  @. prob.u0 = [.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1)]
  prob
end

# set parameters
ρ = 1
γ = 1.1
r = 1.25
e = .2
f = .3
c = 2
d = 2

p = (ρ,γ,r,e,f,c,d)

u0 = [.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1);.1*rand(1)]

genLV = ODEProblem(glv!,u0,tspan,p)
sol = solve(genLV; adaptive=false, dt=dt_tra)

# glv_ensemble = EnsembleProblem(genLV; prob_func = prob_func_glv)
# sim = solve(roe_ensemble; saveat = dt_tra, trajectories = 1000)
tra = Array(sol)

tr = regularize(Dataset(Array(sol)' .+ .0001*randn(size(Array(sol)'))))

figure()
plot(tr[:,1])
plot(tr[:,2])
plot(tr[:,3])
plot(tr[:,4])
plot(tr[:,5])
plot(tr[:,6])
plot(tr[:,7])
plot(tr[:,8])
plot(tr[:,9])


ts = 8

f_slow = 20
f_fast = 8*f_slow
P_fast = 1/f_fast
sampling = P_fast/32
P_slow = 1/f_slow

t = 0:0.0001:1
t = Array(t)
tr = sin.(2*π*t/P_slow).+sin.(2*π*t/P_fast)
tr = regularize(tr .+ .0001*randn(size(tr)))

tr = readdlm("./scripts/generalized Lotka Volterra/data-torus.txt")
tr = Dataset(tr[:,2])
tr = regularize(tr)

ts = 1
w1 = estimate_delay(tr[:,ts], "ac_zero")
# w2 = estimate_delay(tr[:,2], "mi_min")
# w3 = estimate_delay(tr[:,2], "mi_min")
# w4 = estimate_delay(tr[:,2], "mi_min")
# w5 = estimate_delay(tr[:,2], "mi_min")
# w6 = estimate_delay(tr[:,2], "mi_min")
# w7 = estimate_delay(tr[:,2], "mi_min")
# w8 = estimate_delay(tr[:,2], "mi_min")
# w9 = estimate_delay(tr[:,2], "mi_min")
w = w1
w = maximum(hcat(w1,w2,w3,w4,w5,w6,w7,w8,w9))
taus = 0:100

## Perform reconstructions and compute according L-value

# L-value of reference
L_ref = uzal_cost(Dataset(tr), Tw = (4*w), w = w, samplesize=1)

#Standard TDE
@time Y_tde, _ = standard_embedding_cao(tr[:,ts]; method = "ac_zero")
@time Y_tde, _ = standard_embedding_hegger(tr[:,ts]; method = "ac_zero")
L_tde = uzal_cost(Dataset(Y_tde), Tw = (4*w), w = w, samplesize=1)

#MDOP embedding
println("Computation time MDOP:")
@time begin
    tw1 = mdop_maximum_delay(tr[:,ts])

    lm1 = DelayEmbeddings.findlocalminima(tw1[2])
    if lm1[1] ≤ 2
        lmm1 = try lm1[2]
        catch
            lm1[1]
        end
    else
        lmm1 = lm1[1]
    end

    tau_max = lmm1
    taus_mdop = 0:tau_max

    Y_mdop, τ_vals_mdop, ts_vals_mdop, FNNs_mdop , βs_mdop = mdop_embedding(tr[:,ts];
                                                        τs = taus_mdop , w = w)
    L_mdop = uzal_cost(Y_mdop, Tw = (4*w), w = w, samplesize=1)
end

#Garcia & Almeida
println("Computation time Garcia & Almeida method:")
@time Y_GA, τ_vals_GA, ts_vals_GA, FNNs_GA , ns_GA = garcia_almeida_embedding(tr[:,ts];
                                                    τs = taus , w = w, T = w)
L_GA = uzal_cost(Y_GA, Tw = (4*w), w = w, samplesize=1)

#Pecuzal
println("Computation time pecuzal method:")
@time Y_pec, τ_vals_pec, ts_vals_pec, Ls_pec , εs_pec = pecuzal_embedding(tr[:,ts];
                                                            τs = taus , w = w)
L_pec = minimum(Ls_pec)

## compute other evaluation statistics

@time mfnn1, mfnn2, mfnn3, mfnn4, f1, f2, f3, f4, RQA_ref, RQA1, RQA2,
                                            RQA3, RQA4, R_ref, R1, R2, R3, R4 =
        perform_recurrence_analysis(Dataset(tr), Dataset(Y_tde), Dataset(Y_mdop),
                    Dataset(Y_GA), Dataset(Y_pec); ε = 0.08, w = w, kNN = 10)

NN = 500
figure()
subplot(231)
RR_ref = grayscale(R_ref[1:NN,1:NN])
imshow(RR_ref, cmap = "binary_r", extent = (1, size(RR_ref)[1], 1, size(RR_ref)[2]))
title("reference")
subplot(232)
RR1 = grayscale(R1[1:NN,1:NN])
imshow(RR1, cmap = "binary_r", extent = (1, size(RR_ref)[1], 1, size(RR_ref)[2]))
title("TDE")
subplot(233)
RR2 = grayscale(R2[1:NN,1:NN])
imshow(RR2, cmap = "binary_r", extent = (1, size(RR_ref)[1], 1, size(RR_ref)[2]))
title("MDOP")
subplot(234)
RR3 = grayscale(R3[1:NN,1:NN])
imshow(RR3, cmap = "binary_r", extent = (1, size(RR_ref)[1], 1, size(RR_ref)[2]))
title("GA")
subplot(235)
RR4 = grayscale(R4[1:NN,1:NN])
imshow(RR4, cmap = "binary_r", extent = (1, size(RR_ref)[1], 1, size(RR_ref)[2]))
title("PECUZAL")
# ## Analysis of ensemble

figure()
plot3D(Y_tde[:,1],Y_tde[:,2],Y_tde[:,3])
plot3D(Y_mdop[:,1],Y_mdop[:,2],Y_mdop[:,3])
plot3D(Y_pec[:,1],Y_pec[:,2],Y_pec[:,3])

#
RQA_ref.RTE
RQA1.RTE
RQA2.RTE
RQA4.RTE
RQA4.RTE

RQA_ref.DET
RQA1.DET
RQA2.DET
RQA4.DET
