using DrWatson
@quickactivate "PECUZAL_Julia"

using Statistics
using DelimitedFiles
using PyPlot
pygui(true)

## compute mean and std of measures computed in `comm_coupled_roessler.jl`

# In the following the indexing for lines i is as follows:

# i = 1: Garcia & Almeida
# i = 2: MDOP
# i = 3: PECUZAL

μs = 0:0.001:0.1

# load variables
dims = readdlm("./scripts/Fooling systems/synchro results multi 2/dims.csv")
Ls = readdlm("./scripts/Fooling systems/synchro results multi 2/Ls.csv")
τs_GA = readdlm("./scripts/Fooling systems/synchro results multi 2/taus_GA.csv")
τs_mdop = readdlm("./scripts/Fooling systems/synchro results multi 2/taus_mdop.csv")
τs_pec = readdlm("./scripts/Fooling systems/synchro results multi 2/taus_pec.csv")
timewindows = readdlm("./scripts/Fooling systems/synchro results multi 2/timewindows.csv")
ts_GA = readdlm("./scripts/Fooling systems/synchro results multi 2/ts_GA.csv")
ts_mdop = readdlm("./scripts/Fooling systems/synchro results multi 2/ts_mdop.csv")
ts_pec = readdlm("./scripts/Fooling systems/synchro results multi 2/ts_pec.csv")

method = ["GA", "MDOP", "PEC"]
figure()
plot(μs, dims[1,:], label=method[1])
plot(μs, dims[2,:], label=method[2])
plot(μs, dims[3,:], label=method[3])
legend()
grid()
title("reconstruction dimensions")

figure()
subplot(131)
plot(μs, Ls[1,:], label=method[1])
legend()
grid()
title("Ls")
subplot(132)
plot(μs, Ls[2,:], label=method[2])
legend()
grid()
title("Ls")
subplot(133)
plot(μs, Ls[3,:], label=method[3])
legend()
grid()
title("Ls")

figure()
subplot(131)
plot(μs, timewindows[1,:], label=method[1])
legend()
grid()
title("timewindows")
subplot(132)
plot(μs, timewindows[2,:], label=method[2])
legend()
grid()
title("timewindows")
subplot(133)
plot(μs, timewindows[3,:], label=method[3])
legend()
grid()
title("timewindows")
