
using DelayEmbeddings, DynamicalSystemsBase, DelimitedFiles

include("../src/pecora_uzal_method.jl")
include("../src/data_analysis_functions.jl")


s = readdlm("./test/timeseries/lorenz_pecora_uni_x.csv")
s = s[1:3000] # input timeseries = x component of lorenz
w = estimate_delay(s, "mi_min")



# ds = Systems.lorenz96(5; F = 8.0)
# tr = trajectory(ds, 1000.0; Ttr = 100.0, dt = 0.05)
# s = tr[:, 1]
s = readdlm("./test/timeseries/lorenz_pecora_uni_x.csv")
s = vec(s[1:4000]) # input timeseries = x component of lorenz
w = estimate_delay(s, "mi_min")
Tmax = 100
K = 14
samplesize = 1
KNN = 3
Tw = 84

@time Y, τ_vals, ts_vals, Ls , εs = pecuzal_embedding(s;
                                    τs = 0:Tmax , w = w, samplesize = samplesize,
                                    K = K, KNN = KNN, Tw = Tw)



data = readdlm("roessler_test_series.csv")
s = data[1:5000,:]
w1 = estimate_delay(s[:,1], "mi_min")
w2 = estimate_delay(s[:,2], "mi_min")

@time Y_pec, τ_vals_pec, ts_vals_pec, Ls_pec , εs_pec = pecuzal_embedding(Dataset(s[:,2]);
            τs = 0:100, w = 30)

# L1 = uzal_cost(regularize(Dataset(s2[:,2])); samplesize = 1, K = 3, metric = Euclidean(),
#                w = 8)

using PyPlot

pygui(true)

figure()
plot(εs_pec[:,1], label="1")
plot(εs_pec[:,2], label="2")
plot(εs_pec[:,3], label="3")
grid()
