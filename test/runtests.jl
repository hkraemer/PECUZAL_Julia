using DrWatson
@quickactivate "PECUZAL_Julia"
using DelayEmbeddings

ti = time()

const diffeq = (atol = 1e-9, rtol = 1e-9, maxiters = typemax(Int))

include("test_pecuzal_embedding.jl")

ti = time() - ti
println("\nTest took total time of:")
println(round(ti, sigdigits=3), " seconds or ", round(ti/60, sigdigits=3), " minutes")
