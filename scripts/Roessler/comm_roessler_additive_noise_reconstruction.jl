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

# Here we look at a chaotic system: The Roessler attractor in the funnel regime

## Integrate system and determine theiler window

# # set time interval for integration
# N = 5000 # number of samples
# dt_tra = 0.02
# tspan = (0.0, (dt_tra*N))
#
# function roessler!(du,u,p,t)
#     x,y,z = u
#     a,b,c = p
#     du[1] = -(y+z) # dx
#     du[2] = x + (a*y) # dy
#     du[3] = b+z*(x-c) # dz
# end
#
# function multiplicative_noise!(du,u,p,t)
#     du[1] = .3 * u[1]
#     du[2] = .3 * u[2]
#     du[3] = du[3]
# end
#
# a = .2
# b = .2
# c = 5.7
# p = [a,b,c]
#
# function lorenz(du,u,p,t)
#   du[1] = 10.0(u[2]-u[1])
#   du[2] = u[1]*(28.0-u[3]) - u[2]
#   du[3] = u[1]*u[2] - (8/3)*u[3]
# end
#
# # random initial condition
# u0 = [rand(1); rand(1); rand(1)]
#
# prob = SDEProblem(roessler!, multiplicative_noise!, u0, tspan, p)
# sol = solve(prob)
#
# tr = zeros(N,3)
# @inbounds for i = 1:N
#     tr[i,1] = sol.u[i][1]
#     tr[i,2] = sol.u[i][2]
#     tr[i,3] = sol.u[i][3]
# end
#
# figure()
# plot3D(tr[:,1],tr[:,2],tr[:,3])

# roe = Systems.roessler(a = a, b = b, c = c)
roe = Systems.roessler()

# set time interval for integration
N = 5000 # number of samples
dt_tra = 0.05
t = 0:dt_tra:(dt_tra*N)

tr = trajectory(roe, (N*dt_tra); dt = dt_tra, Ttr = (2000*dt_tra))
tr = regularize(tr)

σ = .1
tra = zeros(size(tr))
tra[:,1] = tr[:,1] .+ σ*randn(size(tra,1))
tra[:,2] = tr[:,2] .+ σ*randn(size(tra,1))
tra[:,3] = tr[:,3] .+ σ*randn(size(tra,1))
tr = Dataset(tra)

w1 = estimate_delay(tr[:,1], "mi_min")
w2 = estimate_delay(tr[:,2], "mi_min")

w = maximum(hcat(w1,w2))
taus = 0:100


## Perform reconstructions and compute according L-value

# L-value of reference
L_ref = uzal_cost(tr, Tw = (4*w), w = w, samplesize=1)

#Standard TDE
@time Y_tde, _ = standard_embedding_cao(tr[:,2])
L_tde = uzal_cost(Y_tde, Tw = (4*w), w = w, samplesize=1)

#MDOP embedding
println("Computation time MDOP:")
@time begin
    tw1 = mdop_maximum_delay(tr[:,1])
    tw2 = mdop_maximum_delay(tr[:,2])
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
    tau_max = maximum(hcat(lmm1, lmm2))
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

@time mfnn1, mfnn2, mfnn3, mfnn4, f1, f2, f3, f4, RQA_ref, RQA1, RQA2,
                                            RQA3, RQA4, R_ref, R1, R2, R3, R4 =
        perform_recurrence_analysis(tr, Dataset(Y_tde), Dataset(Y_mdop),
                    Dataset(Y_GA), Dataset(Y_pec); ε = 0.08, w = w, kNN = 10)

NN = 100
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

# figure()
# plot3D(tr[:,1], tr[:,2], tr[:,3])
# title("reference")
# figure()
# plot3D(Y_tde[:,1], Y_tde[:,2], Y_tde[:,3])
# title("TDE")
# figure()
# plot3D(Y_GA[:,1], Y_GA[:,2], Y_GA[:,3])
# title("GA")
# figure()
# plot3D(Y_mdop[:,1], Y_mdop[:,2], Y_mdop[:,3])
# title("MDOP")
# figure()
# plot3D(Y_pec[:,1], Y_pec[:,2], Y_pec[:,3])
# title("PECUZAL")

display("mfnn_tde: $mfnn1")
display("mfnn_mdop: $mfnn2")
display("mfnn_GA: $mfnn3")
display("mfnn_pec: $mfnn4")

display("*****")

display("f_tde: $f1")
display("f_mdop: $f2")
display("f_GA: $f3")
display("f_pec: $f4")

display("*****")

display("DET_ref: $(abs(RQA_ref.DET - RQA_ref.DET)/RQA_ref.DET)")
display("DET_tde: $(abs(RQA_ref.DET - RQA1.DET)/RQA_ref.DET)")
display("DET_mdop: $(abs(RQA_ref.DET - RQA2.DET)/RQA_ref.DET)")
display("DET_GA: $(abs(RQA_ref.DET - RQA3.DET)/RQA_ref.DET)")
display("DET_pec: $(abs(RQA_ref.DET - RQA4.DET)/RQA_ref.DET)")

display("*****")

display("RTE_ref: $(abs(RQA_ref.RTE - RQA_ref.RTE)/RQA_ref.RTE)")
display("RTE_tde: $(abs(RQA_ref.RTE - RQA1.RTE)/RQA_ref.RTE)")
display("RTE_mdop: $(abs(RQA_ref.RTE - RQA2.RTE)/RQA_ref.RTE)")
display("RTE_GA: $(abs(RQA_ref.RTE - RQA3.RTE)/RQA_ref.RTE)")
display("RTE_pec: $(abs(RQA_ref.RTE - RQA4.RTE)/RQA_ref.RTE)")

display("*****")

display("Lmax_ref: $(abs(RQA_ref.Lmax - RQA_ref.Lmax)/RQA_ref.Lmax)")
display("Lmax_tde: $(abs(RQA_ref.Lmax - RQA1.Lmax)/RQA_ref.Lmax)")
display("Lmax_mdop: $(abs(RQA_ref.Lmax - RQA2.Lmax)/RQA_ref.Lmax)")
display("Lmax_GA: $(abs(RQA_ref.Lmax - RQA3.Lmax)/RQA_ref.Lmax)")
display("Lmax_pec: $(abs(RQA_ref.Lmax - RQA4.Lmax)/RQA_ref.Lmax)")

display("*****")

display("ENTR_ref: $(abs(RQA_ref.ENTR - RQA_ref.ENTR)/RQA_ref.ENTR)")
display("ENTR_tde: $(abs(RQA_ref.ENTR - RQA1.ENTR)/RQA_ref.ENTR))")
display("ENTR_mdop: $(abs(RQA_ref.ENTR - RQA2.ENTR)/RQA_ref.ENTR))")
display("ENTR_GA: $(abs(RQA_ref.ENTR - RQA3.ENTR)/RQA_ref.ENTR))")
display("ENTR_pec: $(abs(RQA_ref.ENTR - RQA4.ENTR)/RQA_ref.ENTR))")

display("*****")

display("L_ref: $(abs(RQA_ref.L - RQA_ref.L)/RQA_ref.L)")
display("L_tde: $(abs(RQA_ref.L - RQA1.L)/RQA_ref.L)")
display("L_mdop: $(abs(RQA_ref.L - RQA2.L)/RQA_ref.L)")
display("L_GA: $(abs(RQA_ref.L - RQA3.L)/RQA_ref.L)")
display("L_pec: $(abs(RQA_ref.L - RQA4.L)/RQA_ref.L)")

display("*****")

display("MRT_ref: $(abs(RQA_ref.MRT - RQA_ref.MRT)/RQA_ref.MRT)")
display("MRT_tde: $(abs(RQA_ref.MRT - RQA1.MRT)RQA_ref.MRT)")
display("MRT_mdop: $(abs(RQA_ref.MRT - RQA2.MRT)RQA_ref.MRT)")
display("MRT_GA: $(abs(RQA_ref.MRT - RQA3.MRT)RQA_ref.MRT)")
display("MRT_pec: $(abs(RQA_ref.MRT - RQA4.MRT)RQA_ref.MRT)")

display("*****")

display("LAM_ref: $(abs(RQA_ref.LAM - RQA_ref.LAM)/RQA_ref.LAM)")
display("LAM_tde: $(abs(RQA_ref.LAM - RQA1.LAM)RQA_ref.LAM)")
display("LAM_mdop: $(abs(RQA_ref.LAM - RQA2.LAM)RQA_ref.LAM)")
display("LAM_GA: $(abs(RQA_ref.LAM - RQA3.LAM)RQA_ref.LAM)")
display("LAM_pec: $(abs(RQA_ref.LAM - RQA4.LAM)RQA_ref.LAM)")

display("*****")

display("RTE-MRT_ref: $(abs(RQA_ref.RTE/RQA_ref.MRT - RQA_ref.RTE/RQA_ref.MRT)/(RQA_ref.RTE/RQA_ref.MRT))")
display("RTE-MRT_tde: $(abs(RQA_ref.RTE/RQA_ref.MRT - RQA1.RTE/RQA1.MRT)/(RQA_ref.RTE/RQA_ref.MRT))")
display("RTE-MRT_mdop: $(abs(RQA_ref.RTE/RQA_ref.MRT - RQA2.RTE/RQA2.MRT)/(RQA_ref.RTE/RQA_ref.MRT))")
display("RTE-MRT_GA: $(abs(RQA_ref.RTE/RQA_ref.MRT - RQA3.RTE/RQA3.MRT)/(RQA_ref.RTE/RQA_ref.MRT))")
display("RTE-MRT_pec: $(abs(RQA_ref.RTE/RQA_ref.MRT - RQA4.RTE/RQA4.MRT)/(RQA_ref.RTE/RQA_ref.MRT))")

display("*****")

# # using DelimitedFiles
# #
# writedlm("Duffing_Y_ref",tr)
# writedlm("Duffing_Y_tde",Y_tde)
# writedlm("Duffing_Y_mdop",Y_mdop)
# writedlm("Duffing_Y_GA",Y_GA)
# writedlm("Duffing_Y_pec",Y_pec)
