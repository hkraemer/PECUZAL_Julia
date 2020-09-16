using DrWatson
@quickactivate "new-embedding-methods"

using DelayEmbeddings
using DynamicalSystemsBase
using Test
using Statistics


include("../src/data_analysis_functions.jl")

println("\nTesting fnn_hegger functionality...")

## Try to reproduce Figure 2 in [^Hegger1999]

N = 20000
data = rand(N)
data = regularize(data)

FNN1 = zeros(15)
FNN2 = zeros(15)
FNN3 = zeros(15)
FNN4 = zeros(15)
FNN5 = zeros(15)
FNN6 = zeros(15)

for r = 1:15
    display("r=$r")
    _, fnns, _ = fnn_uniform_hegger(data, 1; r = r, max_dimension = 1, fnn_thres = 0)
    FNN1[r] = fnns[end]
    display(length(fnns))
    _, fnns, _ = fnn_uniform_hegger(data, 1; r = r, max_dimension = 2, fnn_thres = 0)
    FNN2[r] = fnns[end]
    _, fnns, _ = fnn_uniform_hegger(data, 1; r = r, max_dimension = 3, fnn_thres = 0)
    FNN3[r] = fnns[end]
    _, fnns, _ = fnn_uniform_hegger(data, 1; r = r, max_dimension = 4, fnn_thres = 0)
    FNN4[r] = fnns[end]
    _, fnns, _ = fnn_uniform_hegger(data, 1; r = r, max_dimension = 5, fnn_thres = 0)
    FNN5[r] = fnns[end]
    _, fnns, _ = fnn_uniform_hegger(data, 1; r = r, max_dimension = 6, fnn_thres = 0)
    FNN6[r] = fnns[end]
    display(length(fnns))
end



# using PyPlot
# pygui(true)
#
# figure()
# plot(1:15, FNN1, label="m = 1")
# plot(1:15, FNN2, label="m = 2")
# plot(1:15, FNN3, label="m = 3")
# plot(1:15, FNN4, label="m = 4")
# plot(1:15, FNN5, label="m = 5")
# plot(1:15, FNN6, label="m = 6")
# legend()
# ylim(0, 1)
# ylabel("fraction FNNs")
# xlabel("s")
# grid()
