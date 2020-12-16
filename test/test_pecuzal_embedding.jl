using DrWatson
@quickactivate "PECUZAL_Julia"

# These tests are only for DelayEmbeddings.jl are not strictly needed for this
# scientific project

using DelayEmbeddings
using DynamicalSystemsBase
using Test
using DelimitedFiles
using Random

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

    @test -0.728 < Ls[1] < -0.727
    @test -0.323 < Ls[2] < -0.322

    @test τ_vals[2] == 18
    @test τ_vals[3] == 9

    @test length(ts_vals) == 3

end


@testset "Multivariate example" begin
# Test of the proposed Pecora-Uzal-embedding-method (multivariate case)

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
    @test -0.936 < Ls[1] < -0.935
    @test -0.356 < Ls[2] < -0.355
    @test -0.133 < Ls[3] < -0.132
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

    @time Y, τ_vals, ts_vals, Ls , ε★ = pecuzal_embedding(dummy_set;
                                        τs = 0:Tmax , w = w, samplesize = samplesize)

    @test size(Y,2) == 1
end

end
