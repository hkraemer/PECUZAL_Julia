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


idx = 5

Fs = 3.5:0.002:5 # parameter spectrum
F = Fs[idx]
N = 8 # number of oscillators

dt = 0.1 # sampling time
total = 5000  # time series length

# Parameters analysis:
ε = 0.05  # recurrence threshold
dmax = 10   # maximum dimension for traditional tde
lmin = 2   # minimum line length for RQA
trials = 80 # trials for MCDTS
taus = 0:100 # possible delays
Tw = 0  # time window for obtaining the L-value

# randomly pick one time series
t_idx = 2
t_idx = [2,4,7]

# init Lorenz96
u0 = [0.590; 0.766; 0.566; 0.460; 0.794; 0.854; 0.200; 0.298]
lo96 = Systems.lorenz96(N, u0; F = 3.5)

set_parameter!(lo96, 1, F)
data = trajectory(lo96, total*dt; dt = dt, Ttr = 2500 * dt)
data_sample = data[:,t_idx]

# PECUZAL
w = estimate_delay(data_sample[:,2], "mi_min")
# theiler = Int(floor(mean(τ_tde)))
theiler = w
𝒟_pecs, τ_pecs, ts_pecs, Ls_pecs , epss = pecuzal_embedding(data_sample; τs = taus , w = theiler, Tw = theiler)
optimal_d_pecs = size(𝒟_pecs,2)
R = RecurrenceMatrix(𝒟_pec, ε; fixedrate = true)
RQA = rqa(R; theiler = theiler, lmin = lmin)
RQA_pec = hcat(RQA...)
L_pecs = minimum(Ls_pecs)
L_pec[idx]
tau_pec[idx]
