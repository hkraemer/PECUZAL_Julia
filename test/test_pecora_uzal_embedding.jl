using DrWatson
@quickactivate "new-embedding-methods"

using DelayEmbeddings
using DynamicalSystemsBase
#using DifferentialEquations
#using Random
using Test
using DelimitedFiles

include("/Users/hkraemer/Documents/Git/new-embedding-methods/new-embedding-methods/src/pecora_uzal_method.jl")

println("\nTesting pecora_uzal_method.jl...")
@testset "Pecora Uzal" begin
@testset "Univariate example" begin
    ## Test of the proposed Pecora-Uzal-embedding-method (univariate case)

    # For comparison reasons using Travis CI we carry out the integration on a UNIX
    # OS and save the resulting time series
    # lo = Systems.lorenz([1.0, 1.0, 50.0])
    # tr = trajectory(lo, 100; dt = 0.01, Ttr = 10)
    # x = tr[:, 1]
    # writedlm("lorenz_pecora_uni.csv", x)

    s = readdlm("lorenz_pecora_uni.csv")
    s = vec(s) # input timeseries = x component of lorenz
    w = estimate_delay(s, "mi_min")
    Tmax = 100
    K = 14
    samplesize = 1
    KNN = 3
    Tw = 56

    @time Y, τ_vals, ts_vals, Ls , ε★ = pecuzal_embedding(s;
                                        τs = 0:Tmax , w = w, samplesize = samplesize,
                                        K = K, KNN = KNN, Tw = Tw)

    @test -2.949 < Ls[1] < -2.94
    @test -3.099 < Ls[2] < -3.09
    @test -3.039 < Ls[3] < -3.038

    @test τ_vals[2] == 18
    @test τ_vals[3] == 9

    @test length(ts_vals) == 3

end


@testset "Univariate example" begin
## Test of the proposed Pecora-Uzal-embedding-method (multivariate case)

    # For comparison reasons using Travis CI we carry out the integration on a UNIX
    # OS and save the resulting time series
    # lo = Systems.lorenz([1.0, 1.0, 50.0])
    # tr = trajectory(lo, 100; dt = 0.01, Ttr = 10)
    # writedlm("lorenz_pecora_multi.csv", tr)

    tr = readdlm("lorenz_pecora_multi.csv")
    tr = Dataset(tr) # input timeseries = x component of lorenz
    w1 = estimate_delay(s, "mi_min")
    Tmax = 100
    K = 14
    samplesize = 1
    KNN = 3
    Tw = 56

end

end

##

include("/Users/hkraemer/Documents/Git/new-embedding-methods/new-embedding-methods/src/pecora_uzal_method.jl")
tr = readdlm("lorenz_pecora_multi.csv")
tr = Dataset(tr) # input timeseries = x component of lorenz
w1 = estimate_delay(tr[:,1], "mi_min")
w2 = estimate_delay(tr[:,2], "mi_min")
w3 = estimate_delay(tr[:,3], "mi_min")
w = w1
Tmax = 100
K = 14
samplesize = 1
KNN = 3
Tw = 56

@time Y, τ_vals, ts_vals, Ls , ε★ = pecuzal_embedding(tr[1:3000,:];
                                    τs = 0:Tmax , w = w, samplesize = samplesize,
                                    K = K, KNN = KNN, Tw = Tw)
