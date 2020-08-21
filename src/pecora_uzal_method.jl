using DrWatson
@quickactivate "new-embedding-methods"

import Peaks
"""
    pecora_uzal_embedding(s; kwargs...) → Y, τ_vals, ts_vals, Ls ,⟨ε★⟩
A unified approach to properly embed a time series or a set of time series
(`Dataset`) based on the ideas of Pecora et al. [^Pecoral2007] and Uzal et al.
[^Uzal2011].

## Keyword arguments

* `τs= 0:50`: Possible delay values `τs` (in sampling time units). For each of
  the `τs`'s the continuity statistic ⟨ε★⟩ gets computed and further processed
  in order to find optimal delays `τᵢ` for each embedding cycle `i` (read
  algorithm description).
* `w::Int = 1`: Theiler window (neighbors in time with index `w` close to the point,
  that are excluded from being true neighbors). `w=0` means to exclude only the
  point itself, and no temporal neighbors.
* `samplesize::Real = 0.1`: determine the fraction of all phase space points
  (=`length(s)`) to be considered (fiducial points v) to average ε★ to produce
  `⟨ε★⟩`.
* `K::Int = 13`: the amount of nearest neighbors in the δ-ball (read algorithm description).
   Must be at least 8 (in order to gurantee a valid statistic). `⟨ε★⟩` is computed
   taking the minimum result over all `k ∈ K` (read algorithm description).
* `KNN::Int = 3`: the amount of nearest neighbors considered, in order to compute
  σ_k^2 (read algorithm description [`uzal_cost`]@ref). If given a vector, minimum
  result over all `knn ∈ KNN` is returned.
* `Tw::Int = 4*w`: the maximal considered time horizon for obtaining σ_k^2 (read
   algorithm description [`uzal_cost`]@ref).
* `metric = Euclidean()`: metric with which to find nearest neigbhors
* `α::Real = 0.05`: The significance level for obtaining the continuity statistic
* `p::Real = 0.5`: The p-parameter for the binomial distribution used for the
  computation of the continuity statistic ⟨ε★⟩.
* `max_num_of_cycles = 50`: The algorithm will stop after that many cycles no matter what.


## Description
The method works iteratively and gradually builds the final embedding vectors
`Y`. Based on the `⟨ε★⟩`-statistic [`pecora`](@ref) the algorithm picks an
optimal delay value `τᵢ` for each embedding cycle i.
For achieving that, we take the inpute time series `s` and compute the continuity
statistic `⟨ε★⟩`, 1. each local maxima in `⟨ε★⟩` is used for constructing a
candidate embedding trajectory `Y_trial` with a delay corresponding to that
specific peak in `⟨ε★⟩`. 2. We then compute the `L`-statistic [`uzal_cost`](@ref)
for `Y_trial`. 3. We pick the peak/`τ`-value, for which `L` is minimal and
construct the actual embedding trajectory `Y_actual` (1.-3. corresponds to an
embedding cycle). 4. We repeat steps 1.-3. with `Y_actual` as input and stop the
algorithm when `L` can not be reduced anymore. `Y_actual` -> `Y`.

In case of multivariate embedding, i.e. when embedding a set of M time series
(`s::Dataset`), in each embedding cycle `⟨ε★⟩` gets computed for all M time series
available. The optimal delay value `τ` in each embedding cycle is chosen
as the peak/`τ`-value for which `L` is minimal under all available peaks and under
all M `⟨ε★⟩`'s. In the first embedding cycle there will be M! different `⟨ε★⟩`'s
to consider, since it is not clear a priori which time series of the input should
consitute the first component of the embedding vector and form `Y_actual`.

The range of considered delay values is determined in `τs` and for the
nearest neighbor search we respect the Theiler window `w`. The final embedding
vector is stored in `Y` (`Dataset`). The chosen delay values for each embedding
cycle are stored in `τ_vals` and the according time series number chosen for the
each delay value in `τ_vals` is stored in `ts_vals`. For univariate embedding
(`s::Vector`) `ts_vals` is a vector of ones of length `τ_vals`, because there is
simply just one time series to choose from. The function also returns the
`L`-statistic `Ls` for each embedding cycle and the continuity statistic `⟨ε★⟩`
as an `Array` of `Vector`s.

[^Pecora2007]: Pecora, L. M., Moniz, L., Nichols, J., & Carroll, T. L. (2007). [A unified approach to attractor reconstruction. Chaos 17(1)](https://doi.org/10.1063/1.2430294).
[^Uzal2011]: Uzal, L. C., Grinblat, G. L., Verdes, P. F. (2011). [Optimal reconstruction of dynamical systems: A noise amplification approach. Physical Review E 84, 016223](https://doi.org/10.1103/PhysRevE.84.016223).
"""
function pecora_uzal_embedding(s::Vector{T}; τs = 0:50 , w::Int = 1,
    samplesize::Real = 1, K::Int = 13, KNN::Int = 3, Tw::Int=4*w,
    metric = Euclidean(), α::Real = 0.05, p::Real = 0.5,
    max_num_of_cycles = 50) where {T<:Real}

    @assert 0 < samplesize ≤ 1 "Please select a valid `samplesize`, which denotes a fraction of considered fiducial points, i.e. `samplesize` ∈ (0 1]"
    @assert all(x -> x ≥ 0, τs)

    s = regularize(s) # especially important for multivariate embedding and comparison
    # define actual phase space trajectory
    Y_act = Dataset(s)

    # set a flag, in order to tell the while loop when to stop. Each loop
    # stands for encountering a new embedding dimension
    flag, counter = true, 1

    # preallocate output variables
    τ_vals = Int64[0]
    ts_vals = Int64[1]
    Ls = Float64[0.0]
    ε★s = Array{T}(undef, length(τs), max_num_of_cycles)

    # loop over increasing embedding dimensions until some break criterion will
    # tell the loop to stop/break
    while flag
        ε★, _ = pecora(s, Tuple(τ_vals), Tuple(ts_vals); delays = τs, w = w,
                    samplesize = samplesize, K = K, metric = metric, α = α,
                    p = p, undersampling = false)
        ε★s[:, counter] = ε★

        # zero-padding of ⟨ε★⟩ in order to also cover τ=0 (important for the multivariate case)
        ε★ = vec([0; ε★])
        max_idx = Peaks.maxima(ε★) # determine local maxima in ⟨ε★⟩
        L_trials = zeros(Float64, length(max_idx))
        for (i,τ_idx) in enumerate(max_idx)
            # create candidate phase space vector for this peak/τ-value
            Y_trial = DelayEmbeddings.hcat_lagged_values(Y_act,s,τs[τ_idx-1])
            # compute L-statistic
            L_trials[i] = uzal_cost(Y_trial; Tw = Tw, K = KNN, w = w,
                    samplesize = samplesize, metric = metric)
        end

        L_min, min_idx = findmin(L_trials)

        push!(τ_vals, τs[max_idx[min_idx]-1])
        push!(ts_vals, 1)
        counter == 1 ? Ls[1] = L_min : push!(Ls, L_min)

        # create phase space vector for this embedding cycle
        Y_act = DelayEmbeddings.hcat_lagged_values(Y_act,s,τ_vals[counter+1])

        # break criterion 1 (L can not be reduced anymore)
        if counter > 1 && Ls[end]>Ls[end-1]
            println("Algorithm stopped due to minimum L-value reached. VALID embedding achieved.")
            flag = false;
        end
        # break criterion 2 (maximum embedding cycles reached)
        if max_num_of_cycles == counter; flag = false; end

        counter += 1
    end

    return Y_act[:,1:end-1], τ_vals[1:end-1], ts_vals[1:end-1], Ls, ε★s[:,1:counter-1]

end
