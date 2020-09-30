using DrWatson
@quickactivate "new-embedding-methods"

using DifferentialEquations
using DynamicalSystems
using DelayEmbeddings
using DelimitedFiles

include("../../src/pecora_uzal_method.jl")
include("../../src/data_analysis_functions.jl")

Y_GA = readdlm("./scripts/Fooling systems/correlated results/Y_GA.csv")
Y_mdop = readdlm("./scripts/Fooling systems/correlated results/Y_mdop.csv")
Y_pec = readdlm("./scripts/Fooling systems/correlated results/Y_pec.csv")
τ_vals_GA = readdlm("./scripts/Fooling systems/correlated results/taus_GA.csv")
τ_vals_mdop = readdlm("./scripts/Fooling systems/correlated results/taus_mdop.csv")
τ_vals_pec = readdlm("./scripts/Fooling systems/correlated results/taus_pec.csv")
ts_vals_GA = readdlm("./scripts/Fooling systems/correlated results/ts_GA.csv")
ts_vals_mdop = readdlm("./scripts/Fooling systems/correlated results/ts_mdop.csv")
ts_vals_pec = readdlm("./scripts/Fooling systems/correlated results/ts_pec.csv")
