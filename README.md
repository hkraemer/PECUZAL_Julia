# new-embedding-methods

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> new-embedding-methods

It is authored by George Datseris, Hauke Kraemer.

To (locally) reproduce this project, do the following:
0. Download this code base. Notice that raw data are typically not included in the
   git repo and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> cd("path/to/this/project")
   julia> using Pkg; Pkg.activate(".")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.
