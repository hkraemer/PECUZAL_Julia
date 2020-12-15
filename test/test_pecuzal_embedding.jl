using DrWatson
@quickactivate "PECUZAL_Julia"

# These tests are only for DelayEmbeddings.jl are not strictly needed for this
# scientific project

using DelayEmbeddings
using DynamicalSystemsBase
using Test
using DelimitedFiles

include("../src/pecuzal_method.jl")

# println("\nTesting pecuzal_method.jl...")
# @testset "PECUZAL" begin
# @testset "Univariate example" begin
    # For comparison reasons using Travis CI we carry out the integration on a UNIX
    # OS and load the resulting time series (see /test/timeseries/produce_timeseries.jl)

s = readdlm("./test/timeseries/lorenz_pecora_uni_x.csv")
s = vec(s) # input timeseries = x component of lorenz
w = estimate_delay(s, "mi_min")
Tmax = 100
K = 14
samplesize = 1
KNN = 3

@time Y, τ_vals, ts_vals, Ls , εs = pecuzal_embedding(s[1:5000];
                                    τs = 0:Tmax , w = w, samplesize = samplesize,
                                    K = K, KNN = KNN)

@test -0.7264 < Ls[1] < -0.7263
@test -0.3237 < Ls[2] < -0.3236

@test τ_vals[2] == 18
@test τ_vals[3] == 9

@test length(ts_vals) == 3

# end
#
#
# @testset "Multivariate example" begin
## Test of the proposed Pecora-Uzal-embedding-method (multivariate case)

    # For comparison reasons using Travis CI we carry out the integration on a UNIX
    # OS and load the resulting time series (see /test/timeseries/produce_timeseries.jl)

include("../src/pecuzal_method.jl")
tr = readdlm("./test/timeseries/lorenz_pecora_multi.csv")
tr = Dataset(tr) # input timeseries = x component of lorenz
w1 = estimate_delay(tr[:,1], "mi_min")
w2 = estimate_delay(tr[:,2], "mi_min")
w3 = estimate_delay(tr[:,3], "mi_min")
w = w1
Tmax = 100
samplesize = 1

@time Y, τ_vals, ts_vals, Ls , ε★ = pecuzal_embedding(tr[1:5000,:];
                                    τs = 0:Tmax , w = w, samplesize = samplesize)

@test length(ts_vals) == 5
@test ts_vals[3] == ts_vals[4] == ts_vals[5] == 1
@test ts_vals[1] == 3
@test ts_vals[2] == 2
@test τ_vals[2] == 0
@test τ_vals[3] == 62
@test τ_vals[4] == 48
@test τ_vals[5] == 0
@test -0.9362 < Ls[1] < -0.9361
@test -0.3552 < Ls[2] < -0.3551
@test -0.1293 < Ls[3] < -0.1292
@test -0.0144 < Ls[4] < -0.0143

# Dummy input
d1 = randn(1000)
d2 = rand(1000)
dummy_set = Dataset(hcat(d1,d2))

w1 = estimate_delay(d1, "mi_min")
w2 = estimate_delay(d2, "mi_min")
w = w1

@time Y, τ_vals, ts_vals, Ls , ε★ = pecuzal_embedding(dummy_set;
                                    τs = 0:Tmax , w = w, samplesize = samplesize)

@test size(Y,2) == 1
# end
#
# end

include("../src/pecuzal_method.jl")
idx = 487

N = 8 # number of oscillators
Fs = 3.5:0.002:5 # parameter spectrum
dt = 0.1 # sampling time
total = 5000  # time series length
F = Fs[idx]
# Parameters analysis:
ε = 0.05  # recurrence threshold
dmax = 10   # maximum dimension for traditional tde
lmin = 2   # minimum line length for RQA
trials = 80 # trials for MCDTS
taus = 0:100 # possible delays
Tw = 0  # time window for obtaining the L-value

# randomly pick one time series
t_idx = 2
t_idx = [2,4,6]

u0 = [0.590; 0.766; 0.566; 0.460; 0.794; 0.854; 0.200; 0.298]
lo96 = Systems.lorenz96(N, u0; F = 3.5)

set_parameter!(lo96, 1, F)
data = trajectory(lo96, total*dt; dt = dt, Ttr = 2500 * dt)
data_sample = data[:,t_idx]

function space_time_separation(x::Vector{P}, Tw::Int, percentile::Real) where {P}
    N = length(x)
    sts = zeros(Tw+1)
    cnt = 1
    for T = 0:Tw
        sts[cnt] = quantile(abs.(view(x, 1:N-T) .- view(x, 1+T:N)), percentile)
        cnt += 1
    end
    return sts
end

sps = zeros(Tw+1,5)
for (i, perc) in enumerate([.5; .7; .9; .95; .99])
    sps[:,i] = space_time_separation(data_sample, 100, perc)
end
figure()
plot(0:100, sps)
grid()

w = estimate_delay(data_sample, "mi_min")


include("../src/pecuzal_method.jl")
theiler = w
@time 𝒟_pecs, τ_pecs, ts_pecs, Ls_pecs , epss = pecuzal_embedding(data_sample; τs = taus , w = theiler)

Y_act = Dataset(data_sample)
Y_trial = DelayEmbeddings.hcat_lagged_values(Y_act, data_sample, 6)
L1 = zeros(taus[end])
L2 = zeros(taus[end])
dist = zeros(taus[end])
for Tw = 1:taus[end]
    L1[Tw] = uzal_cost(Y_act; Tw = Tw, w = w, samplesize = 1)
    L2[Tw] = uzal_cost(Y_trial; Tw = Tw, w = w, samplesize = 1)
    dist[Tw] = L2[Tw] - L1[Tw]
end

using PyPlot
pygui(true)

figure()
plot(1:taus[end], L1, label="L1")
plot(1:taus[end], L2, label="L2")
plot(1:taus[end], dist, label="Distance")
legend()
grid()


trr = regularize(tr[1:5000,:])
Y_act = DelayEmbeddings.hcat_lagged_values(Dataset(trr[:,2]), trr[:,3], 0)
Y_trial = DelayEmbeddings.hcat_lagged_values(Y_act, trr[:,1], 0)


L1 = zeros(100)
L2 = zeros(100)
for Tw = 1:100
    L1[Tw] = uzal_cost(Y_act; Tw = Tw, w = w, samplesize = 1)
    L2[Tw] = uzal_cost(Y_trial; Tw = Tw, w = w, samplesize = 1)
end

figure()
plot(1:100, L1, label="L1")
plot(1:100, L2, label="L2")
grid()
legend()


figure()
plot(ε★[2][:,1])
grid()


LL1, LL2 = get_minimum_L_by_separation(Y_act, Y_trial, 0:100;
                    w = w)


t = 1:1000
data = sin.(2*π*t/60)
w = estimate_delay(data, "mi_min")
theiler = w
𝒟_pecs, τ_pecs, ts_pecs, Ls_pecs , epss = pecuzal_embedding(data; τs = 0:50 , w = theiler)
