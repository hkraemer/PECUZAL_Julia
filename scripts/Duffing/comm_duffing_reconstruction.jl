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

# Here we look at a limit cycle example: the Duffing oscillator

## Integrate system and determine theiler window

duf = Systems.duffing(ω = 1, f = 0.3, d = 0.2, β = -1)
# set time interval for integration
N = 10000 # number of samples
dt_tra = 0.5
t = 0:dt_tra:(dt_tra*N)


function duffing!(du,u,p,t)
    x,y,z = u
    δ,α,β,γ,ω = p

    du[1] = y # dx
    du[2] = -δ*y - α*x - β*x^3 + z #dy
    du[3] = γ*ω*cos(ω*t)# dz
end

# set parameters for chaotic Duffing system
δ = .2
α = 1
β = -1
γ = .3
ω = 1
p = [δ,α,β,γ,ω]

# random initial condition
u0 = [2*rand(1); 2*rand(1); 2*rand(1)]
# construct a ContinuousDynamicalSystem
duffing_system = ContinuousDynamicalSystem(duffing!, u0, p)
# solve ODEproblem / get trajectory
tr = trajectory(duffing_system, (N*dt_tra);  Ttr = (2000*dt_tra))


tr = trajectory(duf, (N*dt_tra); dt = dt_tra, Ttr = (2000*dt_tra))
tr = regularize(tr)

w1 = estimate_delay(tr[:,1], "mi_min")
w2 = estimate_delay(tr[:,2], "mi_min")
w3 = estimate_delay(tr[:,3], "mi_min")

w = maximum(hcat(w1,w2,w3))
taus = 0:200

## Perform reconstructions and compute according L-value

# L-value of reference
L_ref = uzal_cost(tr, Tw = (4*w), w = w, samplesize=1)

#Standard TDE
@time Y_tde, _ = standard_embedding_cao(tr[:,1])
L_tde = uzal_cost(Y_tde, Tw = (4*w), w = w, samplesize=1)


#MDOP embedding
println("Computation time MDOP:")
@time begin
    tw1 = mdop_maximum_delay(tr[:,1])
    tw2 = mdop_maximum_delay(tr[:,2])
    #tw3 = mdop_maximum_delay(tr[:,3])
    #tau_max = maximum(hcat(tw1[1], tw2[1], tw3[1]))
    tau_max = maximum(hcat(tw1[1], tw2[1]))
    taus_mdop = 0:tau_max

    Y_mdop, τ_vals_mdop, ts_vals_mdop, FNNs_mdop , βs_mdop = mdop_embedding(tr[:,1:2];
                                                        τs = taus_mdop , w = w)
    L_mdop = uzal_cost(Y_mdop, Tw = (4*w), w = w, samplesize=1)
end

#Garcia & Almeida
println("Computation time Garcia & Almeida method:")
@time Y_GA, τ_vals_GA, ts_vals_GA, FNNs_GA , ns_GA = garcia_almeida_embedding(tr[:,1:2];
                                                    τs = taus , w = w, T = w)
L_GA = uzal_cost(Y_GA, Tw = (4*w), w = w, samplesize=1)

#Pecuzal
println("Computation time pecuzal method:")
@time Y_pec, τ_vals_pec, ts_vals_pec, Ls_pec , εs_pec = pecuzal_embedding(tr[:,1:2];
                                                            τs = taus , w = w)
L_pec = minimum(Ls_pec)

## compute other evaluation statistics

@time f1, f2, f3, f4, RQA_ref, RQA1, RQA2, RQA3, RQA4 = perform_recurrence_analysis(
                        tr, Y_tde, Y_mdop, Y_GA, Y_pec;
                        ε = 0.08, w_ref = w, w₁ = w,  w₂ = w, w₃ = w, w₄ = w)

display("f_tde: $f1")
display("f_mdop: $f2")
display("f_GA: $f3")
display("f_pec: $f4")

display("DET_ref: $(abs(RQA_ref.DET - RQA_ref.DET))")
display("DET_tde: $(abs(RQA_ref.DET - RQA1.DET))")
display("DET_mdop: $(abs(RQA_ref.DET - RQA2.DET))")
display("DET_GA: $(abs(RQA_ref.DET - RQA3.DET))")
display("DET_pec: $(abs(RQA_ref.DET - RQA4.DET))")

display("RTE_ref: $(abs(RQA_ref.RTE - RQA_ref.RTE))")
display("RTE_tde: $(abs(RQA_ref.RTE - RQA1.RTE))")
display("RTE_mdop: $(abs(RQA_ref.RTE - RQA2.RTE))")
display("RTE_GA: $(abs(RQA_ref.RTE - RQA3.RTE))")
display("RTE_pec: $(abs(RQA_ref.RTE - RQA4.RTE))")

# display("Lmax_ref: $(abs(RQA_ref.Lmax - RQA_ref.Lmax))")
# display("Lmax_tde: $(abs(RQA_ref.Lmax - RQA1.Lmax))")
# display("Lmax_mdop: $(abs(RQA_ref.Lmax - RQA2.Lmax))")
# display("Lmax_GA: $(abs(RQA_ref.Lmax - RQA3.Lmax))")
# display("Lmax_pec: $(abs(RQA_ref.Lmax - RQA4.Lmax))")
#
# display("ENTR_ref: $(abs(RQA_ref.ENTR - RQA_ref.ENTR))")
# display("ENTR_tde: $(abs(RQA_ref.ENTR - RQA1.ENTR))")
# display("ENTR_mdop: $(abs(RQA_ref.ENTR - RQA2.ENTR))")
# display("ENTR_GA: $(abs(RQA_ref.ENTR - RQA3.ENTR))")
# display("ENTR_pec: $(abs(RQA_ref.ENTR - RQA4.ENTR))")
#
# display("L_ref: $(abs(RQA_ref.L - RQA_ref.L))")
# display("L_tde: $(abs(RQA_ref.L - RQA1.L))")
# display("L_mdop: $(abs(RQA_ref.L - RQA2.L))")
# display("L_GA: $(abs(RQA_ref.L - RQA3.L))")
# display("L_pec: $(abs(RQA_ref.L - RQA4.L))")
#
# display("MRT_ref: $(abs(RQA_ref.MRT - RQA_ref.MRT))")
# display("MRT_tde: $(abs(RQA_ref.MRT - RQA1.MRT))")
# display("MRT_mdop: $(abs(RQA_ref.MRT - RQA2.MRT))")
# display("MRT_GA: $(abs(RQA_ref.MRT - RQA3.MRT))")
# display("MRT_pec: $(abs(RQA_ref.MRT - RQA4.MRT))")

display("RTE-MRT_ref: $(abs(RQA_ref.RTE/RQA_ref.MRT - RQA_ref.RTE/RQA_ref.MRT))")
display("RTE-MRT_tde: $(abs(RQA_ref.RTE/RQA_ref.MRT - RQA1.RTE/RQA1.MRT))")
display("RTE-MRT_mdop: $(abs(RQA_ref.RTE/RQA_ref.MRT - RQA2.RTE/RQA2.MRT))")
display("RTE-MRT_GA: $(abs(RQA_ref.RTE/RQA_ref.MRT - RQA3.RTE/RQA3.MRT))")
display("RTE-MRT_pec: $(abs(RQA_ref.RTE/RQA_ref.MRT - RQA4.RTE/RQA4.MRT))")


using PyPlot

pygui(true)

figure()
plot(tr[:,1], tr[:,2])
