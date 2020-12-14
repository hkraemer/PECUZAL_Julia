using DrWatson
@quickactivate "PECUZAL_Julia"

# These tests are only for DelayEmbeddings.jl are not strictly needed for this
# scientific project

using DelayEmbeddings
using DynamicalSystemsBase
using Test
using DelimitedFiles

include("../src/pecuzal_method.jl")

println("\nTesting pecuzal_method.jl...")
@testset "PECUZAL" begin
@testset "Univariate example" begin
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

    @test -3.049 < Ls[1] < -3.048
    @test -3.366 < Ls[2] < -3.365
    @test -3.3201 < Ls[3] < -3.3200

    @test τ_vals[2] == 18
    @test τ_vals[3] == 9

    @test length(ts_vals) == 3

end


@testset "Multivariate example" begin
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


    @test length(ts_vals) == 3
    @test ts_vals[1] == 2
    @test ts_vals[2] == 3
    @test ts_vals[3] == 1
    @test τ_vals[1] == 0
    @test τ_vals[2] == 0
    @test τ_vals[3] == 0
    @test -2.984 < Ls[1] < -2.983
    @test -3.264 < Ls[2] < -3.263
    @test -3.153 < Ls[3] < -3.152

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
end

end


##
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

theiler = w
𝒟_pecs, τ_pecs, ts_pecs, Ls_pecs , epss = pecuzal_embedding(data_sample; τs = taus , w = theiler, Tw = 5)

include("../src/pecuzal_method.jl")


Tw = 100
E2 = zeros(8, Tw)
epss = zeros(8)
for (i, fid) in enumerate(sample(1:length(data_sample), 8; replace = false))
    println(fid)
    _, E2[i,:], epss[i] = uzal_cost2(Dataset(data_sample); w = theiler, samplesize=1, Tw=Tw)
end
means = zeros(8,Tw)
medians = zeros(8,Tw)
for i = 1:Tw
    means[:,i] = mean(E2[:,1:i], dims=2)
    medians[:,i] = median(E2[:,1:i], dims=2)
end

figure(figsize=(20,10))
for i = 1:8
    subplot(2,4,i)
    plot(1:Tw, E2[i,:], label="E₂")
    plot(1:Tw, means[i,:], linestyle = "dotted", label="mean")
    plot(1:Tw, medians[i,:], linestyle = "dashed", label="median")
    xlabel("Tw")
    ylabel("E₂ value")
    legend()
    grid()
end

figure(figsize=(20,10))
for i = 1:8
    subplot(2,4,i)
    plot(1:Tw, means[i,:]/epss[i], label="mean")
    plot(1:Tw, medians[i,:]/epss[i], linestyle = "dotted",label="median")
    title("normalized to neighbourhood")
    xlabel("Tw")
    ylabel("E₂ value")
    legend()
    grid()
end

using Plots

maxis , maxi_idx = get_maxima(vec(epss))

YY = DelayEmbeddings.hcat_lagged_values(Dataset(data_sample), data_sample, 6)
Lss = zeros(100)
Lss2 = zeros(100)
for Tw = 1:100
    Lss[Tw] = mean(DelayEmbeddings.uzal_cost_local(Dataset(data_sample); w = theiler, samplesize=1, Tw=Tw))
    Lss2[Tw] = mean(DelayEmbeddings.uzal_cost_local(YY; w = theiler, samplesize=1, Tw=Tw))
    # Lss[Tw] = DelayEmbeddings.uzal_cost(Dataset(data_sample); w = theiler, samplesize=1, Tw=Tw)
    # Lss2[Tw] = DelayEmbeddings.uzal_cost(YY; w = theiler, samplesize=1, Tw=Tw)
end

using Plots
plot(Lss, label="1")
plot!(Lss2)


plot(epss)
xgrid!
ygrid!

Tw = 34
Lss = DelayEmbeddings.uzal_cost_local(Dataset(data_sample); w = theiler, samplesize=1, Tw=Tw)
Lss2 = DelayEmbeddings.uzal_cost_local(YY; w = theiler, samplesize=1, Tw=Tw)
histogram(Lss)
histogram!(Lss2)
title!("2-dim Tw = $Tw")
