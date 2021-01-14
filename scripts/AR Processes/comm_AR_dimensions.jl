using DrWatson
@quickactivate "PECUZAL_Julia"

using DelayEmbeddings
using DelimitedFiles
using Random

include("../../src/pecuzal_method.jl")
include("../../src/data_analysis_functions.jl")

# compute time series
Random.seed!(1234)

u0_1 = randn(3)
AR3 = ar_process3(u0_1, .45, 0.3, .1, 1.0, 40000)
#AR3 = ar_process3(u0_1, .45, 0.1, .4, 1.0, 40000)

u0_2 = randn(5)
AR5 = ar_process5(u0_2, .2, 0.1, 0.05, .3, 0.2, 1.0, 40000)

u0_3 = randn(7)
AR7 = ar_process7(u0_3, .2, 0.04, 0.0, 0.05, .1, 0.2, .3, 1.0, 40000)
#AR7 = ar_process7(u0_3, .1, 0.04, 0.0, 0.05, .3, 0.2, .3, 1.0, 40000)

figure()
subplot(311)
plot(AR3)
subplot(312)
plot(AR5)
subplot(313)
plot(AR7)

ts_lengths = [500, 1000, 3000, 6000, 10000, 15000, 20000, 30000, 40000]

dim_AR3 = zeros(length(ts_lengths))
dim_AR5 = zeros(length(ts_lengths))
dim_AR7 = zeros(length(ts_lengths))
Ls_3 = []
Ls_5 = []
Ls_7 = []
τ_vals_3 = []
τ_vals_5 = []
τ_vals_7 = []

for (idx,ts_length) in enumerate(ts_lengths)
    println("TS-length: $ts_length")
    w3 = estimate_delay(AR3[1:ts_length],"mi_min")
    w5 = estimate_delay(AR5[1:ts_length],"mi_min")
    w7 = estimate_delay(AR7[1:ts_length],"mi_min")

    Y3, τ3, _, Ls_3_, _ = pecuzal_embedding_update(AR3[1:ts_length]; τs = 0:50 , w = w3)
    Y5, τ5, _, Ls_5_, _ = pecuzal_embedding_update(AR5[1:ts_length]; τs = 0:50 , w = w3)
    Y7, τ7, _, Ls_7_, _ = pecuzal_embedding_update(AR7[1:ts_length]; τs = 0:50 , w = w7)

    push!(Ls_3, Ls_3_)
    push!(Ls_5, Ls_5_)
    push!(Ls_7, Ls_7_)
    push!(τ_vals_3, τ3)
    push!(τ_vals_5, τ5)
    push!(τ_vals_7, τ7)
    dim_AR3[idx] = size(Y3,2)
    dim_AR5[idx] = size(Y5,2)
    dim_AR7[idx] = size(Y7,2)
    if idx == 5
        break
    end
end

# writedlm("./scripts/AR Processes/AR3_dims.csv", dim_AR3)
# writedlm("./scripts/AR Processes/AR5_dims.csv", dim_AR5)
# writedlm("./scripts/AR Processes/AR7_dims.csv", dim_AR7)
