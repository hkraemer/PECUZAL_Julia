# using DrWatson
# @quickactivate "new-embedding-methods"

using DelayEmbeddings
using DynamicalSystemsBase
using Random

export pecuzal_embedding

"""
    pecuzal_embedding(s; kwargs...) → Y, τ_vals, ts_vals, Ls ,⟨ε★⟩
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
function pecuzal_embedding(s::Vector{T}; τs = 0:50 , w::Int = 1,
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
    Ls = Float64[]
    ε★s = Array{T}(undef, length(τs), max_num_of_cycles)

    # loop over increasing embedding dimensions until some break criterion will
    # tell the loop to stop/break
    while flag
        Y_act = pecuzal_embedding_cycle!(
                Y_act, flag, s, τs, w, counter, ε★s, τ_vals, metric,
                Ls, ts_vals, samplesize, K, α, p, Tw, KNN)

        flag = pecuzal_break_criterion(Ls, counter, max_num_of_cycles)
        counter += 1
    end

    return Y_act[:,1:end-1], τ_vals[1:end-1], ts_vals[1:end-1], Ls, ε★s[:,1:counter-1]

end

function pecuzal_embedding(Y::Dataset{D, T}; τs = 0:50 , w::Int = 1,
    samplesize::Real = 1, K::Int = 13, KNN::Int = 3, Tw::Int=4*w,
    metric = Euclidean(), α::Real = 0.05, p::Real = 0.5,
    max_num_of_cycles = 50) where {D, T<:Real}

    @assert 0 < samplesize ≤ 1 "Please select a valid `samplesize`, which denotes a fraction of considered fiducial points, i.e. `samplesize` ∈ (0 1]"
    @assert all(x -> x ≥ 0, τs)

    # centralize data
    Y = regularize(Y)

    # define actual phase space trajectory
    Y_act = []

    # set a flag, in order to tell the while loop when to stop. Each loop
    # stands for encountering a new embedding dimension
    flag, counter = true, 1

    # preallocate output variables
    τ_vals = Int64[0]
    ts_vals = Int64[]
    Ls = Float64[]
    #ε★s = Array{T}(undef, length(τs), max_num_of_cycles)
    ε★s = fill(zeros(T, length(τs), size(Y,2)), 1, max_num_of_cycles)

    # loop over increasing embedding dimensions until some break criterion will
    # tell the loop to stop/break
    while flag
        Y_act = pecuzal_multivariate_embedding_cycle!(
                Y_act, flag, Y, τs, w, counter, ε★s, τ_vals, metric,
                Ls, ts_vals, samplesize, K, α, p, Tw, KNN)

        flag = pecuzal_break_criterion(Ls, counter, max_num_of_cycles)
        counter += 1
    end

    return Y_act[:,1:end-1], τ_vals[1:end-1], ts_vals[1:end-1], Ls, ε★s[:,1:counter-1]

end


"""
Perform one embedding cycle on `Y`
"""
function pecuzal_embedding_cycle!(
        Y, flag, s, τs, w, counter, ε★s, τ_vals, metric,
        Ls, ts_vals, samplesize, K, α, p, Tw, KNN)

    ε★, _ = pecora(s, Tuple(τ_vals), Tuple(ts_vals); delays = τs, w = w,
                samplesize = samplesize, K = K, metric = metric, α = α,
                p = p, undersampling = false)
    ε★s[:,counter] = ε★

    # zero-padding of ⟨ε★⟩ in order to also cover τ=0 (important for the multivariate case)
    ε★ = vec([0; ε★])
    # get the L-statistic for each peak in ⟨ε★⟩ and take the one according to L_min
    L_trials, max_idx = local_L_statistics(ε★, Y, s, τs, Tw, KNN, w, samplesize, metric)
    L_min, min_idx = findmin(L_trials)

    push!(τ_vals, τs[max_idx[min_idx]-1])
    push!(ts_vals, 1)
    push!(Ls, L_min)

    # create phase space vector for this embedding cycle
    Y_act = DelayEmbeddings.hcat_lagged_values(Y,s,τ_vals[counter+1])

    return Y_act
end

"""
Perform one embedding cycle on `Y` with a multivariate set Ys
"""
function pecuzal_multivariate_embedding_cycle!(
        Y, flag, Ys, τs, w, counter, ε★s, τ_vals, metric,
        Ls, ts_vals, samplesize, K, α, p, Tw, KNN)

    M = size(Ys,2)
    # in the 1st cycle we have to check all (size(Y,2))! combinations
    if counter == 1
        L_min = zeros(M)
        L_min_idx = zeros(Int, M)
        ε★ = zeros(length(τs), M*M)
        idx = zeros(Int, M)
        for ts = 1:M
            ε★[:,1+(M*(ts-1)):M*ts], _ = pecora(Ys, (0,), (ts,); delays = τs,
                        w = w, samplesize = samplesize, K = K, metric = metric,
                        α = α, p = p, undersampling = false)
            L_min[ts], L_min_idx[ts], idx[ts] = choose_right_embedding_params(
                                            ε★[:,1+(M*(ts-1)):M*ts], Ys[:,ts],
                                            Ys, τs, Tw, KNN, w, samplesize,
                                            metric)
        end
        L_mini, min_idx = findmin(L_min)
        # update τ_vals, ts_vals, Ls, ε★s
        push!(τ_vals, τs[L_min_idx[min_idx]])
        push!(ts_vals, min_idx)             # time series to start with
        push!(ts_vals, idx[min_idx])        # result of 1st embedding cycle
        push!(Ls, L_mini)
        ε★s[counter] = ε★[:,1+(M*(ts_vals[1]-1)):M*ts_vals[1]]

        # create phase space vector for this embedding cycle
        Y_act = DelayEmbeddings.hcat_lagged_values(Ys[:,ts_vals[counter]],
                                     Ys[:,ts_vals[counter+1]],τ_vals[counter+1])

    # in all other cycles we just have to check (size(Y,2)) combinations
    else
        ε★, _ = pecora(Ys, Tuple(τ_vals), Tuple(ts_vals); delays = τs, w = w,
                samplesize = samplesize, K = K, metric = metric, α = α,
                p = p, undersampling = false)
        # update τ_vals, ts_vals, Ls, ε★s
        choose_right_embedding_params!(ε★, Y, Ys, τ_vals, ts_vals, Ls, ε★s,
                                    counter, τs, Tw, KNN, w, samplesize, metric)
        # create phase space vector for this embedding cycle
        Y_act = DelayEmbeddings.hcat_lagged_values(Y, Ys[:, ts_vals[counter+1]],
                                                            τ_vals[counter+1])
    end

    return Y_act
end

"""
Choose the minimum L and corresponding τ for each ε★-statistic
"""
function choose_right_embedding_params!(ε★, Y, Ys, τ_vals, ts_vals, Ls, ε★s,
                                 counter, τs, Tw, KNN, w, samplesize, metric)
    L_min_ = zeros(size(Ys,2))
    τ_idx = zeros(Int,size(Ys,2))
    for ts = 1:size(Ys,2)
        # zero-padding of ⟨ε★⟩ in order to also cover τ=0 (important for the multivariate case)
        # get the L-statistic for each peak in ⟨ε★⟩ and take the one according to L_min
        L_trials_, max_idx_ = local_L_statistics(vec([0; ε★[:,ts]]), Y, Ys[:,ts],
                                        τs, Tw, KNN, w, samplesize, metric)
        L_min_[ts], min_idx_ = findmin(L_trials_)
        τ_idx[ts] = max_idx_[min_idx_]-1
    end
    idx = sortperm(L_min_)
    L_mini, min_idx = findmin(L_min_)
    push!(τ_vals, τs[τ_idx[min_idx]])
    push!(ts_vals, min_idx)
    push!(Ls, L_mini)

    ε★s[counter] = ε★
end

function choose_right_embedding_params(ε★, Y, Ys, τs, Tw, KNN, w, samplesize, metric)
    L_min_ = zeros(size(Ys,2))
    τ_idx = zeros(Int,size(Ys,2))
    for ts = 1:size(Ys,2)
        # zero-padding of ⟨ε★⟩ in order to also cover τ=0 (important for the multivariate case)
        # get the L-statistic for each peak in ⟨ε★⟩ and take the one according to L_min
        L_trials_, max_idx_ = local_L_statistics(vec([0; ε★[:,ts]]), Y, Ys[:,ts],
                                        τs, Tw, KNN, w, samplesize, metric)
        L_min_[ts], min_idx_ = findmin(L_trials_)
        τ_idx[ts] = max_idx_[min_idx_]-1
    end
    idx = sortperm(L_min_)
    return L_min_[idx[1]], τ_idx[idx[1]], idx[1]

end


"""
Return the L-statistic `L` and indices `max_idx` for all local maxima in ε★
"""
function local_L_statistics(ε★, Y, s, τs, Tw, KNN, w, samplesize, metric)
    max_idx = get_maxima(ε★) # determine local maxima in ⟨ε★⟩
    L_trials = zeros(Float64, length(max_idx))
    for (i,τ_idx) in enumerate(max_idx)
        # create candidate phase space vector for this peak/τ-value
        Y_trial = DelayEmbeddings.hcat_lagged_values(Y,s,τs[τ_idx-1])
        # compute L-statistic
        L_trials[i] = uzal_cost(Y_trial; Tw = Tw, K = KNN, w = w,
                samplesize = samplesize, metric = metric)
    end
    return L_trials, max_idx
end

function pecuzal_break_criterion(Ls, counter, max_num_of_cycles)
    flag = true
    if counter > 1 && Ls[end]>Ls[end-1]
        println("Algorithm stopped due to minimum L-value reached. "*
                "VALID embedding achieved ✓.")
        flag = false
    end
    if max_num_of_cycles == counter
        println("Algorithm stopped due to hitting max cycle number. "*
                "Valid embedding NOT achieved ⨉.")
        flag = false
    end
    return flag
end


"""
Return the maxima of the given time series s
"""
function get_maxima(s::Vector{T}) where {T}
    maximas = Int[]
    N = length(s)
    flag = false
    first_point = 0
    for i = 2:N-1
        if s[i-1] < s[i] && s[i+1] < s[i]
            flag = false
            push!(maximas, i)
        end
        # handling constant values
        if flag
            if s[i+1] < s[first_point]
                flag = false
                push!(maximas, first_point)
            elseif s[i+1] > s[first_point]
                flag = false
            end
        end
        if s[i-1] < s[i] && s[i+1] == s[i]
            flag = true
            first_point = i
        end
    end
    # make sure there is no empty vector returned
    if isempty(maximas)
        _, maximas = findmax(s)
    end
    return maximas
end
