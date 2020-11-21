using DrWatson
@quickactivate "PECUZAL_Julia"
using DelayEmbeddings

# These tests are only for DelayEmbeddings.jl are not strictly needed for this
# scientific project

ti = time()

include("test_pecuzal_embedding.jl")

ti = time() - ti
println("\nTest took total time of:")
println(round(ti, sigdigits=3), " seconds or ", round(ti/60, sigdigits=3), " minutes")
