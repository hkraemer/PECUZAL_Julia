using DelayEmbeddings
using DynamicalSystemsBase
using Test
using DelimitedFiles
using PyPlot
pygui(true)

# using DrWatson
# @quickactivate "new-embedding-methods"

include("../src/pecora_uzal_method.jl")

N = 300

tr = readdlm("./test/timeseries/lorenz_pecora_multi.csv")
tr = Dataset(tr) # input timeseries = x component of lorenz
tra = tr[1:N,:]
w1 = estimate_delay(tr[:,1], "mi_min")
w2 = estimate_delay(tr[:,2], "mi_min")
w3 = estimate_delay(tr[:,3], "mi_min")
w = w1
w = 10
Tmax = 100
samplesize = 1
deltas = 14

## start debug
taus = 0:Tmax
@time Y, τ_vals, ts_vals, Ls , ε★ = pecuzal_embedding(tr[1:N,:];
                                    τs = 0:Tmax , w = w, samplesize = samplesize,
                                    K = deltas)

epsis1 = pecora(tra, (0,), (1,); delays = 0:Tmax, samplesize = samplesize, w = w,
                                            metric = Euclidean(), K = deltas)
epsis2 = pecora(tra, (0,), (2,); delays = 0:Tmax, samplesize = samplesize, w = w,
                                        metric = Euclidean(), K = deltas)
epsis3 = pecora(tra, (0,), (3,); delays = 0:Tmax, samplesize = samplesize, w = w,
                                            metric = Euclidean(), K = deltas)

## by hand!
delays = 0:Tmax
tra = regularize(tra)
vspace = genembed(tra, (0,), (1,))
vtree = KDTree(vspace.data[1:end-maximum(delays)], Euclidean())
samplesize=1
L = length(vspace)
if samplesize==1
    ns = vec(1:length(vec(max(1, (-minimum(delays) + 1)):min(L, L - maximum(delays)))))
else
    ns = sample(vec(max(1, (-minimum(delays) + 1)):min(L, L - maximum(delays))),
    length(vec(max(1, (-minimum(delays) + 1)):min(L, L - maximum(delays)))),
    replace = false)
end
vs = vspace[ns]
allNNidxs, allNNdist = DelayEmbeddings.all_neighbors(vtree, vs, ns, 20, 1)

x = vec(Matrix(regularize(vspace)))
epsis_raw_1 = DelayEmbeddings.continuity_per_timeseries(x, ns, allNNidxs, delays, 13, 0.05, 0.5)


display("from x: $(epsis1[1][1:4,1])")
display("from y: $(epsis2[1][1:4,1])")
display("from z: $(epsis3[1][1:4,1])")

figure()
subplot(221)
plot(taus, ε★[1])
grid()
title("pecuzal 1st")
subplot(222)
plot(taus, epsis1[1])
grid()
title("pec x")
subplot(223)
plot(taus, epsis2[1])
grid()
title("pec y")
subplot(224)
plot(taus, epsis3[1])
grid()
title("pec z")






## end debug


@time Y, τ_vals, ts_vals, Ls , ε★ = pecuzal_embedding(tr[1:5000,:];
                                    τs = 0:Tmax , w = w, samplesize = samplesize)


@test length(ts_vals) == 4
@test ts_vals[3] == ts_vals[4] == 1
@test ts_vals[1] == 2
@test ts_vals[2] == 3
@test τ_vals[2] == 0
@test τ_vals[3] == 62
@test τ_vals[4] == 0
@test -2.19 < Ls[1] < -2.18
@test -2.49 < Ls[2] < -2.48
@test -2.57 < Ls[3] < -2.56

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

Y = regularize(Y)
L = uzal_cost(Y; Tw = 4*w, w=w, samplesize=samplesize)


using PyPlot

pygui(true)

figure()
plot(0:Tmax,ε★[1])
grid()
