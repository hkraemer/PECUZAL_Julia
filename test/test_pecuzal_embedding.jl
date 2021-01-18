using DrWatson
@quickactivate "PECUZAL_Julia"

# These tests are only for development and further merging the functionality into
# DelayEmbeddings.jl and are not strictly needed for this scientific project

using DelayEmbeddings
using DynamicalSystemsBase
using Test
using DelimitedFiles
using Random

include("../src/pecuzal_method.jl")
include("../src/data_analysis_functions.jl")

println("\nTesting pecuzal_method.jl for comparison to MATLAB and Python implementations:")
@testset "PECUZAL other algorithm comparison" begin
    s = readdlm("./test/timeseries/lorenz_pecora_uni_x.csv")
    s = vec(s[1:500]) # input timeseries = x component of lorenz
    Tmax = 100
    theiler = 21

    @time Y, τ_vals, ts_vals, Ls , εs = pecuzal_embedding_update(s;
                               τs = 0:Tmax , w = theiler)
    @test size(Y,2) == 4
    @test -0.6134 < sum(Ls) < -0.6133
    @test τ_vals[2]  == 21
    @test τ_vals[3]  == 13
    @test τ_vals[4]  == 78

    @time Y, τ_vals, ts_vals, Ls , εs = pecuzal_embedding_update(s;
                               τs = 0:Tmax , w = theiler, econ = true)

    @test -0.6078 < sum(Ls) < -0.6077
    @test τ_vals[2]  == 21
    @test τ_vals[3]  == 13
    @test τ_vals[4]  == 78

    tr = readdlm("./test/timeseries/lorenz_pecora_multi.csv")
    tr = Dataset(tr) # input timeseries = x component of lorenz
    w = 15
    Tmax = 100

    @time Y, τ_vals, ts_vals, Ls , ε★ = pecuzal_embedding_update(tr[1:500,1:2];
                                         τs = 0:Tmax , w = w)
    @test size(Y,2) == 2
    @test ts_vals[1] == 2
    @test ts_vals[2] == 1
    @test τ_vals[1] == τ_vals[2] == 0
    @test sum(Ls) < -0.5505736

    @time Y, τ_vals, ts_vals, Ls , ε★ = pecuzal_embedding_update(tr[1:500,1:2];
                                         τs = 0:Tmax , w = w, econ = true)
    @test size(Y,2) == 2
    @test ts_vals[1] == 2
    @test ts_vals[2] == 1
    @test τ_vals[1] == τ_vals[2] == 0
    @test sum(Ls) < -0.544942
end

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

    @time Y, τ_vals, ts_vals, Ls , εs = pecuzal_embedding_update(s[1:5000];
                                        τs = 0:Tmax , w = w, samplesize = samplesize,
                                        K = K, KNN = KNN)
    L1 = sum(Ls)
    L2 = compute_delta_L(s[1:5000], τ_vals, Tmax;
         KNN = KNN, w = w)
    @test L1 == L2
    @test -0.728 < Ls[1] < -0.727
    @test -0.3234 < Ls[2] < -0.3233

    @test τ_vals[2] == 18
    @test τ_vals[3] == 9

    @test length(ts_vals) == 3

    @time Y, τ_vals, ts_vals, Ls , εs = pecuzal_embedding_update(s[1:5000];
                                        τs = 0:Tmax , w = w, samplesize = samplesize,
                                        K = K, KNN = KNN, econ = true)

    L1 = sum(Ls)

    @test L1 < -1.047465
    @test -0.7263 < Ls[1] < -0.7261
    @test -0.3213 < Ls[2] < -0.3212

    @test τ_vals[2] == 18
    @test τ_vals[3] == 9

    @test length(ts_vals) == 3

    YY1 = DelayEmbeddings.hcat_lagged_values(Dataset(s), s, 21)
    YY2 = DelayEmbeddings.hcat_lagged_values(YY1, s, 13)
    @test length(YY1) == length(YY2)
    @test YY1 == YY2[:,1:2]

end


@testset "Multivariate example" begin
# Test of the proposed Pecora-Uzal-embedding-method (multivariate case)

    # For comparison reasons using Travis CI we carry out the integration on a UNIX
    # OS and load the resulting time series (see /test/timeseries/produce_timeseries.jl)

    include("../src/pecuzal_method.jl")
    include("../src/data_analysis_functions.jl")
    tr = readdlm("./test/timeseries/lorenz_pecora_multi.csv")
    tr = Dataset(tr) # input timeseries = x component of lorenz
    w1 = estimate_delay(tr[:,1], "mi_min")
    w2 = estimate_delay(tr[:,2], "mi_min")
    w3 = estimate_delay(tr[:,3], "mi_min")
    w = w1
    Tmax = 100
    samplesize = 1
    KNN = 3

    @time Y, τ_vals, ts_vals, Ls , ε★ = pecuzal_embedding_update(tr[1:5000,:];
                                         τs = 0:Tmax , w = w, samplesize = samplesize)
    L1 = sum(Ls)
    L2 = compute_delta_L(tr[1:5000,:], τ_vals, ts_vals, Tmax; KNN = KNN, w = w)
    @test L1 == L2
    @test length(ts_vals) == 5
    @test ts_vals[3] == ts_vals[4] == ts_vals[5] == 1
    @test ts_vals[1] == 3
    @test ts_vals[2] == 2
    @test τ_vals[2] == 0
    @test τ_vals[3] == 62
    @test τ_vals[4] == 48
    @test τ_vals[5] == 0
    @test -0.936 < Ls[1] < -0.935
    @test -0.356 < Ls[2] < -0.355
    @test -0.131 < Ls[3] < -0.130
    @test -0.015 < Ls[4] < -0.014

    # Dummy input
    include("../src/pecuzal_method.jl")
    Random.seed!(1234)
    d1 = randn(1000)
    d2 = rand(1000)
    dummy_set = Dataset(hcat(d1,d2))
    w1 = estimate_delay(d1, "mi_min")
    w2 = estimate_delay(d2, "mi_min")
    w = w1

    @time Y, τ_vals, ts_vals, Ls , ε★ = pecuzal_embedding_update(dummy_set;
                                        τs = 0:Tmax , w = w, samplesize = samplesize)

    @test size(Y,2) == 1

end

end
