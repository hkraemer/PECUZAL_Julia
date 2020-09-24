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

# Here we try to fool the reconstruction algorithms by feeding in time series
# from the Roessler system and the Lorenz system as well as internally
# correlated values of these time series.

## Generate data

# set time interval for integration
N = 5000 # number of samples
dt_tra = 0.03
transients = 2000
tspan = (dt_tra*N) + (dt_tra*transients)

function coupled_roessler!(du,u,p,t)
    x1,x2,x3,y1,y2,y3 = u
    a,b,c,ν,μ = p
    du[1] = -(1 + ν)*x2 - x3
    du[2] = (1 + ν)*x1 + a*x2 + μ*(y2-x2)
    du[3] = b + x3*(x1-c)
    du[4] = -(1 - ν)*y2 - y3
    du[5] = (1 - ν)*y1 + a*y2 + μ*(x2-y2)
    du[6] = b + y3*(y1-c)
end


a = .16
b = .1
c = 8.5
ν = 0.02

μs = 0:0.001:0.1
Ls = zeros(3,length(μs))
timewindows = zeros(3,length(μs))
dims = zeros(3,length(μs))
ws = zeros(length(μs))
ts_pec = []
ts_GA = []
ts_mdop = []
τs_pec = []
τs_GA = []
τs_mdop = []

taus = 0:100

u0 = [2*rand(1); 2*rand(1); 2*rand(1); 2*rand(1); 2*rand(1); 2*rand(1)]
@time for (cnt,μ) in enumerate(μs)
    p = [a,b,c,ν,μ]
    roe = ODEProblem(coupled_roessler!,u0,tspan,p)
    sol = solve(roe; saveat = dt_tra)
    tr = Array(sol)'[transients+1:end,:]
    tr = regularize(Dataset(tr))

    w1 = estimate_delay(tr[:,1], "mi_min")
    w2 = estimate_delay(tr[:,2], "mi_min")
    w3 = estimate_delay(tr[:,3], "mi_min")
    w4 = estimate_delay(tr[:,4], "mi_min")
    w5 = estimate_delay(tr[:,5], "mi_min")
    w6 = estimate_delay(tr[:,6], "mi_min")

    w = maximum(hcat(w1,w2,w3,w4,w5,w6))

    ws[cnt] = w

    if cnt == 1
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
        global taus_mdop = 0:tau_max
    end

    #Garcia & Almeida
    Y_GA, τ_vals_GA, ts_vals_GA, FNNs_GA , ns_GA = garcia_almeida_embedding(tr;
                                                        τs = taus , w = w, T = w)
    timewindows[1,cnt] = sum(τ_vals_GA)
    dims[1,cnt] = size(Y_GA,2)
    Ls[1,cnt] = uzal_cost(Y_GA, Tw = (4*w), w = w, samplesize=1)
    push!(ts_GA, ts_vals_GA)
    push!(τs_GA, τ_vals_GA)

    #MDOP
    Y_mdop, τ_vals_mdop, ts_vals_mdop, FNNs_mdop , βs_mdop = mdop_embedding(tr;
                                                        τs = taus_mdop , w = w)
    timewindows[2,cnt] = sum(τ_vals_mdop)
    dims[2,cnt] = size(Y_mdop,2)
    Ls[2,cnt] = uzal_cost(Y_mdop, Tw = (4*w), w = w, samplesize=1)
    push!(ts_mdop, ts_vals_mdop)
    push!(τs_mdop, τ_vals_mdop)

    #Pecuzal
    Y_pec, τ_vals_pec, ts_vals_pec, Ls_pec , εs_pec = pecuzal_embedding(tr;
                                                                τs = taus , w = w)
    timewindows[3,cnt] = sum(τ_vals_pec)
    dims[3,cnt] = size(Y_pec,2)
    Ls[3,cnt] = minimum(Ls_pec)
    push!(ts_pec, ts_vals_pec)
    push!(τs_pec, τ_vals_pec)
end

writedlm("./scripts/Fooling systems/synchro results multi 1/Ls.csv", Ls)
writedlm("./scripts/Fooling systems/synchro results multi 1/timewindows.csv", timewindows)
writedlm("./scripts/Fooling systems/synchro results multi 1/dims.csv", dims)
writedlm("./scripts/Fooling systems/synchro results multi 1/ws.csv", ws)
writedlm("./scripts/Fooling systems/synchro results multi 1/ts_pec.csv", ts_pec)
writedlm("./scripts/Fooling systems/synchro results multi 1/ts_GA.csv", ts_GA)
writedlm("./scripts/Fooling systems/synchro results multi 1/ts_mdop.csv", ts_mdop)
writedlm("./scripts/Fooling systems/synchro results multi 1/taus_pec.csv", τs_pec)
writedlm("./scripts/Fooling systems/synchro results multi 1/taus_GA.csv", τs_GA)
writedlm("./scripts/Fooling systems/synchro results multi 1/taus_mdop.csv", τs_mdop)


## multi 2

@time for (cnt,μ) in enumerate(μs)
    p = [a,b,c,ν,μ]
    roe = ODEProblem(coupled_roessler!,u0,tspan,p)
    sol = solve(roe; saveat = dt_tra)
    tr = Array(sol)'[transients+1:end,:]
    tr = regularize(Dataset(tr))

    w2 = estimate_delay(tr[:,2], "mi_min")
    w5 = estimate_delay(tr[:,5], "mi_min")

    w = maximum(hcat(w2,w5))
    ws[cnt] = w

    tra = zeros(size(tr,1),2)
    tra[:,1] = tr[:,2]
    tra[:,2] = tr[:,2]
    tra = Dataset(tra)


    if cnt == 1
        tw1 = mdop_maximum_delay(tra[:,1])
        tw2 = mdop_maximum_delay(tra[:,2])

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
        global taus_mdop = 0:tau_max
    end

    #Garcia & Almeida
    Y_GA, τ_vals_GA, ts_vals_GA, FNNs_GA , ns_GA = garcia_almeida_embedding(tra;
                                                        τs = taus , w = w, T = w)
    timewindows[1,cnt] = sum(τ_vals_GA)
    dims[1,cnt] = size(Y_GA,2)
    Ls[1,cnt] = uzal_cost(Y_GA, Tw = (4*w), w = w, samplesize=1)
    push!(ts_GA, ts_vals_GA)
    push!(τs_GA, τ_vals_GA)

    #MDOP
    Y_mdop, τ_vals_mdop, ts_vals_mdop, FNNs_mdop , βs_mdop = mdop_embedding(tra;
                                                        τs = taus_mdop , w = w)
    timewindows[2,cnt] = sum(τ_vals_mdop)
    dims[2,cnt] = size(Y_mdop,2)
    Ls[2,cnt] = uzal_cost(Y_mdop, Tw = (4*w), w = w, samplesize=1)
    push!(ts_mdop, ts_vals_mdop)
    push!(τs_mdop, τ_vals_mdop)

    #Pecuzal
    Y_pec, τ_vals_pec, ts_vals_pec, Ls_pec , εs_pec = pecuzal_embedding(tra;
                                                                τs = taus , w = w)
    timewindows[3,cnt] = sum(τ_vals_pec)
    dims[3,cnt] = size(Y_pec,2)
    Ls[3,cnt] = minimum(Ls_pec)
    push!(ts_pec, ts_vals_pec)
    push!(τs_pec, τ_vals_pec)
end

writedlm("./scripts/Fooling systems/synchro results multi 2/Ls.csv", Ls)
writedlm("./scripts/Fooling systems/synchro results multi 2/timewindows.csv", timewindows)
writedlm("./scripts/Fooling systems/synchro results multi 2/dims.csv", dims)
writedlm("./scripts/Fooling systems/synchro results multi 2/ws.csv", ws)
writedlm("./scripts/Fooling systems/synchro results multi 2/ts_pec.csv", ts_pec)
writedlm("./scripts/Fooling systems/synchro results multi 2/ts_GA.csv", ts_GA)
writedlm("./scripts/Fooling systems/synchro results multi 2/ts_mdop.csv", ts_mdop)
writedlm("./scripts/Fooling systems/synchro results multi 2/taus_pec.csv", τs_pec)
writedlm("./scripts/Fooling systems/synchro results multi 2/taus_GA.csv", τs_GA)
writedlm("./scripts/Fooling systems/synchro results multi 2/taus_mdop.csv", τs_mdop)
