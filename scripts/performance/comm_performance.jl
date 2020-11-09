using DrWatson
@quickactivate "new-embedding-methods"

using DynamicalSystems
using DelayEmbeddings
using DelimitedFiles
using BenchmarkTools

include("../../src/pecuzal_method.jl")
include("../../src/data_analysis_functions.jl")

## We analyze the computational complexity of the proposed PECUZAL method in
# comparison to TDE, the method from G&A and MDOP. We feed in the y-component
# of the Rössler system for different time series length and compute the average
# computation time.

roe = Systems.roessler()

# set time interval for integration
N = 20000 # number of samples
dt_tra = 0.05
t = 0:dt_tra:(dt_tra*N)

tr = trajectory(roe, (N*dt_tra); dt = dt_tra, Ttr = (2000*dt_tra))
tr = regularize(tr)

# define function which wrap the whole embedding process
function tde_time(data, w, taus)
    Y_tde, τ_tde = standard_embedding_cao(data; τs = taus)
end

function GA_time(data, w, taus)
    Y_GA, τ_vals_GA, ts_vals_GA, FNNs_GA , ns_GA = garcia_almeida_embedding(data;
                                                        τs = taus , w = w, T = w)
end

function mdop_time(data, w)
    tw1 = mdop_maximum_delay(data)
    lm1 = DelayEmbeddings.findlocalminima(tw1[2])
    if lm1[1] ≤ 2
        lmm1 = try lm1[2]
        catch
            lm1[1]
        end
    else
        lmm1 = lm1[1]
    end

    taus_mdop = 0:lmm1

    Y_mdop, τ_vals_mdop, ts_vals_mdop, FNNs_mdop , βs_mdop = mdop_embedding(data;
                                                        τs = taus_mdop , w = w)
end

function pec_time(data, w, taus)
    Y_pec, τ_vals_pec, ts_vals_pec, Ls_pec , εs_pec = pecuzal_embedding(data;
                                                                τs = taus , w = w)
end

step = 1000
Nss = 1000:step:N

times_tde = zeros(length(Nss))
times_GA = zeros(length(Nss))
times_mdop = zeros(length(Nss))
times_pec = zeros(length(Nss))

for (i,Ns) in enumerate(Nss)
    display("run: $i")
    data = tr[1:Ns,2]
    # compute Theiler window
    w = estimate_delay(data, "mi_min", 0:150)
    # and range of encountered time delays
    taus1 = 1:(4*w)
    taus2 = 0:(4*w)

    b_tde = @benchmark tde_time($data, $w, $taus1)
    b_GA = @benchmark GA_time($data, $w, $taus2)
    b_mdop = @benchmark mdop_time($data, $w)
    b_pec = @benchmark pec_time($data, $w, $taus2)

    times_tde[i] = median(b_tde).time
    times_GA[i] = median(b_GA).time
    times_mdop[i] = median(b_mdop).time
    times_pec[i] = median(b_pec).time

end

writedlm("./scripts/performance/results/times_tde.csv", times_tde)
writedlm("./scripts/performance/results/times_GA.csv", times_GA)
writedlm("./scripts/performance/results/times_mdop.csv", times_mdop)
writedlm("./scripts/performance/results/times_pec.csv", times_pec)
writedlm("./scripts/performance/results/times.csv", Nss)
