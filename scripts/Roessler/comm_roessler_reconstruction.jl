using DelayEmbeddings

#Download some test timeseries
repo = "https://raw.githubusercontent.com/JuliaDynamics/ExercisesRepo/master/timeseries"
tsfolder = joinpath(@__DIR__, "timeseries")
todownload = ["$n.csv" for n in 1:4]

mkpath(tsfolder)
for a in todownload
    download(repo*"/"*a, joinpath(tsfolder, a))
end

ti = time()

const diffeq = (atol = 1e-9, rtol = 1e-9, maxiters = typemax(Int))
