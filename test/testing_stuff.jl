using DelayEmbeddings
using DynamicalSystemsBase
using Test
using DelimitedFiles


include("../src/pecora_uzal_method.jl")

# For comparison reasons using Travis CI we carry out the integration on a UNIX
# OS and load the resulting time series (see /test/timeseries/produce_timeseries.jl)

s = readdlm("./test/timeseries/lorenz_pecora_uni.csv")
s = vec(s) # input timeseries = x component of lorenz
w = estimate_delay(s, "mi_min")
Tmax = 100
K = 14
samplesize = 1
KNN = 3
Tw = 56

@time Y, τ_vals, ts_vals, Ls , εs = pecuzal_embedding(s[1:5000];
                                    τs = 0:Tmax , w = w, samplesize = samplesize,
                                    K = K, KNN = KNN, Tw = Tw)

@test -2.73 < Ls[1] < -2.72
@test -2.77 < Ls[2] < -2.76
@test -2.73 < Ls[3] < -2.72

@test τ_vals[2] == 18
@test τ_vals[3] == 9

@test length(ts_vals) == 3
