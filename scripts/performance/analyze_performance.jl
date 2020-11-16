using DrWatson
@quickactivate "PECUZAL_Julia"

using DelimitedFiles
using PyPlot
pygui(true)

times_tde = readdlm("./scripts/performance/results/times_tde.csv")
times_GA = readdlm("./scripts/performance/results/times_GA.csv")
times_mdop = readdlm("./scripts/performance/results/times_mdop.csv")
times_pec = readdlm("./scripts/performance/results/times_pec.csv")
Nss = readdlm("./scripts/performance/results/times.csv")

## Plot results

lw = 2.5      # linewidth
ms = 7        # markersize
fsa = 14      # fontsize
axislabelsize = 10
lwa = 1.5                 # axis line width
ticklabelsize = 10
fsl = 10        # legendfontsize

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
font0 = Dict(
        "font.size" => fsa,
        "font.weight" => "bold",
        "axes.labelweight" => "bold",
        "axes.titleweight" => "bold",
        "axes.labelsize" => axislabelsize,
        "axes.linewidth" => lwa,
        "xtick.labelsize" => ticklabelsize,
        "ytick.labelsize" => ticklabelsize,
        "legend.fontsize" => fsl)
merge!(rcParams, font0)

figure()
plot(Nss, times_tde./1000000, label="TDE", marker="o", linestyle="solid", linewidth=lw, markersize=ms)
plot(Nss, times_GA/1000000, label="G&A", marker="D", linestyle="solid", linewidth=lw, markersize=ms)
plot(Nss, times_mdop/1000000, label="MDOP", marker="s", linestyle="solid", linewidth=lw, markersize=ms)
plot(Nss, times_pec/1000000, label="PECUZAL", marker="x", linestyle="solid", linewidth=lw, markersize=ms)
plot(Nss, Nss.*log.(Nss), label="N log(N)", linestyle="dashed")
legend()
xlabel("time series length N")
ylabel("time [μs]")
title("  Median computation time \n of reconstruction methods")
yscale("log")
grid()
