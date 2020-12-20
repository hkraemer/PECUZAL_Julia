using DrWatson
@quickactivate "PECUZAL_Julia"

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

# Here we try to fool the reconstruction algorithms by feeding in time series
# from the Roessler system and the Lorenz system as well as internally
# correlated values of these time series.

## Generate data

# set time interval for integration
N = 5000 # number of samples
dt_tra = 0.03
transients = 2000
tspan = (dt_tra*N) + (dt_tra*transients)

roe = Systems.roessler()
lo = Systems.lorenz()

tr1 = trajectory(roe, (N*dt_tra); dt = dt_tra, Ttr = (transients*dt_tra))
tr2 = trajectory(lo, (N*dt_tra); dt = dt_tra, Ttr = (transients*dt_tra))

tr = zeros(N+1,6)
tr[:,1] = tr1[:,1]
tr[:,2] = tr1[:,2]
tr[:,3] = tr2[:,1]
tr[:,4] = tr2[:,2]
tr[:,5] = (tr1[:,1].*tr1[:,1]).^2 .+ tr1[:,2]
tr[:,6] = (tr2[:,2].*tr2[:,1]).^2 .+ tr2[:,2]

tr = regularize(Dataset(tr))

method = "mi_min"

w1 = estimate_delay(tr[:,1], method)
w2 = estimate_delay(tr[:,2], method)
w3 = estimate_delay(tr[:,3], method)
w4 = estimate_delay(tr[:,4], method)
w5 = estimate_delay(tr[:,5], method)
w6 = estimate_delay(tr[:,6], method)

w = maximum(hcat(w1,w2,w3,w4,w5,w6))
taus = 0:100

# L-value of reference
L_ref = uzal_cost(tr, Tw = (4*w), w = w, samplesize=1)

#MDOP embedding
println("Computation time MDOP:")
@time begin
    tw1 = mdop_maximum_delay(tr[:,1])
    tw2 = mdop_maximum_delay(tr[:,2])
    tw3 = mdop_maximum_delay(tr[:,3])
    tw4 = mdop_maximum_delay(tr[:,4])
    tw5 = mdop_maximum_delay(tr[:,5])
    tw6 = mdop_maximum_delay(tr[:,6])
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
    lm4 = DelayEmbeddings.findlocalminima(tw4[2])
    if lm4[1] ≤ 2
        lmm4 = try lm4[2]
        catch
            lm4[1]
        end
    else
        lmm4 = lm4[1]
    end
    lm5 = DelayEmbeddings.findlocalminima(tw5[2])
    if lm5[1] ≤ 2
        lmm5 = try lm5[2]
        catch
            lm5[1]
        end
    else
        lmm5 = lm5[1]
    end
    lm6 = DelayEmbeddings.findlocalminima(tw6[2])
    if lm6[1] ≤ 2
        lmm6 = try lm6[2]
        catch
            lm6[1]
        end
    else
        lmm6 = lm6[1]
    end
    tau_max = maximum(hcat(lmm1, lmm2, lmm3, lmm4, lmm5, lmm6))
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
@time Y_pec, τ_vals_pec, ts_vals_pec, Ls_pec , εs_pec = pecuzal_embedding_update(tr;
                                                            τs = taus , w = w)
L_pec = minimum(Ls_pec)


using PyPlot
pygui(true)

figure()
plot(Y_GA[:,1], Y_GA[:,2], Y_GA[:,3])
title("GA")

figure()
plot3D(Y_mdop[:,1], Y_mdop[:,2], Y_mdop[:,3])
title("mdop")

figure()
plot3D(Y_pec[:,1], Y_pec[:,2], Y_pec[:,3])
title("pec")

# writedlm("./scripts/Fooling systems/correlated results/Y_GA.csv",Y_GA)
# writedlm("./scripts/Fooling systems/correlated results/Y_mdop.csv",Y_mdop)
# writedlm("./scripts/Fooling systems/correlated results/Y_pec.csv",Y_pec)
# writedlm("./scripts/Fooling systems/correlated results/taus_GA.csv", τ_vals_GA)
# writedlm("./scripts/Fooling systems/correlated results/taus_mdop.csv", τ_vals_mdop)
# writedlm("./scripts/Fooling systems/correlated results/taus_pec.csv", τ_vals_pec)
# writedlm("./scripts/Fooling systems/correlated results/ts_GA.csv", ts_vals_GA)
# writedlm("./scripts/Fooling systems/correlated results/ts_mdop.csv", ts_vals_mdop)
# writedlm("./scripts/Fooling systems/correlated results/ts_pec.csv", ts_vals_pec)
