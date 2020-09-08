## This script produces the time series we use in the tests

using DifferentialEquations
using Random
using DelimitedFiles

# For comparison reasons using Travis CI we carry out the integration on a UNIX
# OS and save the resulting time series

# Lorenz
lo = Systems.lorenz([1.0, 1.0, 50.0])
tr = trajectory(lo, 100; dt = 0.01, Ttr = 10)
x = tr[:, 1] # x-component of time series
y = tr[:, 2] #y-component of time series
writedlm("./test/timeseries/lorenz_pecora_uni_x.csv", x)
writedlm("./test/timeseries/lorenz_pecora_uni_y.csv", y)
writedlm("./test/timeseries/lorenz_pecora_multi.csv", tr)


# Lorenz2
lo = Systems.lorenz([1.0, 1.0, 50.0]; ρ=60)
tr = trajectory(lo, 100; dt = 0.01, Ttr = 10)
x = tr[:, 1] # x-component of time series
y = tr[:, 2] #y-component of time series
writedlm("./test/timeseries/lorenz2_pecora_uni_x.csv", x)
writedlm("./test/timeseries/lorenz2_pecora_uni_y.csv", y)
writedlm("./test/timeseries/lorenz2_pecora_multi.csv", tr)
