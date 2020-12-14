using DrWatson
@quickactivate "PECUZAL_Julia"

# Core functions for the usage in `/scripts`-folder
using DelayEmbeddings
using DynamicalSystemsBase
using Random
using Neighborhood
using StatsBase
using Distances

export pecuzal_embedding

"""
    pecuzal_embedding(s; kwargs...) → Y, τ_vals, ts_vals, Ls ,⟨ε★⟩
A unified approach to properly embed a time series or a set of time series
(`Dataset`) based on the ideas of Pecora et al. [^Pecoral2007] and Uzal et al.
[^Uzal2011].

## Keyword arguments

* `τs = 0:50`: Possible delay values `τs` (in sampling time units). For each of
  the `τs`'s the continuity statistic ⟨ε★⟩ gets computed and further processed
  in order to find optimal delays `τᵢ` for each embedding cycle `i` (read
  algorithm description).
* `w::Int = 1`: Theiler window (neighbors in time with index `w` close to the point,
  that are excluded from being true neighbors). `w=0` means to exclude only the
  point itself, and no temporal neighbors.
* `samplesize::Real = 1`: determine the fraction of all phase space points
  (=`length(s)`) to be considered (fiducial points v) to average ε★, in order to
  produce `⟨ε★⟩`.
* `K::Int = 13`: the amount of nearest neighbors in the δ-ball (read algorithm
  description). Must be at least 8 (in order to gurantee a valid statistic).
  `⟨ε★⟩` is computed taking the minimum result over all `k ∈ K`.
* `KNN::Int = 3`: the amount of nearest neighbors considered, in order to compute
  σ_k^2 (read algorithm description [`uzal_cost`]@ref). If given a vector, the
  minimum result over all `knn ∈ KNN` is returned.
* `α::Real = 0.05`: The significance level for obtaining the continuity statistic
* `p::Real = 0.5`: The p-parameter for the binomial distribution used for the
  computation of the continuity statistic ⟨ε★⟩.
* `max_cycles = 50`: The algorithm will stop after that many cycles no matter what.


## Description
The method works iteratively and gradually builds the final embedding vectors
`Y`. Based on the `⟨ε★⟩`-statistic [`pecora`](@ref) the algorithm picks an
optimal delay value `τᵢ` for each embedding cycle `i`.
For achieving that, we take the inpute time series `s`, denoted as the actual
phase space trajectory `Y_actual` and compute the continuity statistic `⟨ε★⟩`.
1. Each local maxima in `⟨ε★⟩` is used for constructing a
candidate embedding trajectory `Y_trial` with a delay corresponding to that
specific peak in `⟨ε★⟩`. 2. We then compute the `L`-statistic [`uzal_cost`](@ref)
for `Y_trial` (`L-trial`) and `Y_actual` (`L_actual`) for increasing prediction time
horizons (free parameter in the `L`-statistic) and save the maximum difference
`max(L-trial - L_actual)` as `L_decrease` (Note that this is a negative number,
since the `L`-statistic decreases with better reconstructions). 3. We pick the
peak/`τ`-value, for which `L_decrease` is minimal (=maximum decrease of the overall
`L`-value) and construct the actual embedding trajectory `Y_actual`
(steps 1.-3. correspond to an embedding cycle). 4. We repeat steps 1.-3. with
`Y_actual` as input and stop the algorithm when `L_decrease` is > 0, i.e. when
and additional embedding component would not lead to a lower overall L-value.
`Y_actual` -> `Y`.

In case of multivariate embedding, i.e. when embedding a set of M time series
(`s::Dataset`), in each embedding cycle the continuity statistic `⟨ε★⟩` gets
computed for all M time series available. The optimal delay value `τ` in each
embedding cycle is chosen as the peak/`τ`-value for which `L_decrease` is
minimal under all available peaks and under all M `⟨ε★⟩`'s. In the first
embedding cycle there will be M^2 different `⟨ε★⟩`'s to consider, since it is
not clear a priori which time series of the input should consitute the first
component of the embedding vector and form `Y_actual`.

The range of considered delay values is determined in `τs` and for the
nearest neighbor search we respect the Theiler window `w`. The final embedding
vector is stored in `Y` (`Dataset`). The chosen delay values for each embedding
cycle are stored in `τ_vals` and the according time series numbers chosen for
each delay value in `τ_vals` are stored in `ts_vals`. For univariate embedding
(`s::Vector`) `ts_vals` is a vector of ones of length `τ_vals`, because there is
simply just one time series to choose from. The function also returns the
`L_decrease`-values for each embedding cycle and the continuity statistic `⟨ε★⟩`
as an `Array` of `Vector`s.

For distance computations the Euclidean norm is used.

[^Pecora2007]: Pecora, L. M., Moniz, L., Nichols, J., & Carroll, T. L. (2007). [A unified approach to attractor reconstruction. Chaos 17(1)](https://doi.org/10.1063/1.2430294).
[^Uzal2011]: Uzal, L. C., Grinblat, G. L., Verdes, P. F. (2011). [Optimal reconstruction of dynamical systems: A noise amplification approach. Physical Review E 84, 016223](https://doi.org/10.1103/PhysRevE.84.016223).
"""
function pecuzal_embedding(s::Vector{T}; τs = 0:50 , w::Int = 1,
    samplesize::Real = 1, K::Int = 13, KNN::Int = 3,
    α::Real = 0.05, p::Real = 0.5, max_cycles::Int = 50) where {T<:Real}

    @assert 0 < samplesize ≤ 1 "Please select a valid `samplesize`, which denotes a fraction of considered fiducial points, i.e. `samplesize` ∈ (0 1]"
    @assert all(x -> x ≥ 0, τs)
    metric = Euclidean()

    s_orig = s
    s = regularize(s) # especially important for comparative L-statistics
    # define actual phase space trajectory
    Y_act = Dataset(s)

    # set a flag, in order to tell the while loop when to stop. Each loop
    # stands for encountering a new embedding dimension
    flag, counter = true, 1

    # preallocate output variables
    τ_vals = Int64[0]
    ts_vals = Int64[1]
    Ls = Float64[]
    ε★s = Array{T}(undef, length(τs), max_cycles)

    # loop over increasing embedding dimensions until some break criterion will
    # tell the loop to stop/break
    while flag
        Y_act = pecuzal_embedding_cycle!(
                Y_act, flag, s, τs, w, counter, ε★s, τ_vals, metric,
                Ls, ts_vals, samplesize, K, α, p, KNN)
        flag = pecuzal_break_criterion(Ls, counter, max_cycles)
        counter += 1
    end
    # construct final reconstruction vector
    NN = (length(s)-sum(τ_vals[1:end-1]))
    Y_final = s_orig
    for i = 2:length(τ_vals[1:end-1])
        Y_final = DelayEmbeddings.hcat_lagged_values(Y_final, s_orig, τ_vals[i])
    end
    return Y_final, τ_vals[1:end-1], ts_vals[1:end-1], Ls[1:end-1], ε★s[:,1:counter-1]

end

function pecuzal_embedding(Y::Dataset{D, T}; τs = 0:50 , w::Int = 1,
    samplesize::Real = 1, K::Int = 13, KNN::Int = 3,
    α::Real = 0.05, p::Real = 0.5, max_cycles::Int = 50) where {D, T<:Real}

    @assert 0 < samplesize ≤ 1 "Please select a valid `samplesize`, which denotes a fraction of considered fiducial points, i.e. `samplesize` ∈ (0 1]"
    @assert all(x -> x ≥ 0, τs)
    metric = Euclidean()

    Y_orig = Y
    Y = regularize(Y) # especially important for comparative L-statistics

    # define actual phase space trajectory
    Y_act = []

    # set a flag, in order to tell the while loop when to stop. Each loop
    # stands for encountering a new embedding dimension
    flag, counter = true, 1

    # preallocate output variables
    τ_vals = Int64[0]
    ts_vals = Int64[]
    Ls = Float64[]
    ε★s = fill(zeros(T, length(τs), size(Y,2)), 1, max_cycles)

    # loop over increasing embedding dimensions until some break criterion will
    # tell the loop to stop/break
    while flag
        Y_act = pecuzal_multivariate_embedding_cycle!(
                Y_act, flag, Y, τs, w, counter, ε★s, τ_vals, metric,
                Ls, ts_vals, samplesize, K, α, p, KNN)

        flag = pecuzal_break_criterion(Ls, counter, max_cycles)
        counter += 1
    end
    # construct final reconstruction vector
    Y_final = Y_orig[:,ts_vals[1]]
    for i = 2:length(τ_vals[1:counter-1])
        Y_final = DelayEmbeddings.hcat_lagged_values(Y_final,Y_orig[:,ts_vals[i]],τ_vals[i])
    end

    return Y_final, τ_vals[1:end-1], ts_vals[1:end-1], Ls[1:end-1], ε★s[:,1:counter-1]

end


"""
Perform one univariate embedding cycle on `Y_act`. Return the new `Y_act`
"""
function pecuzal_embedding_cycle!(
        Y_act, flag, s, τs, w, counter, ε★s, τ_vals, metric,
        Ls, ts_vals, samplesize, K, α, p, KNN)

    ε★, _ = pecora(s, Tuple(τ_vals), Tuple(ts_vals); delays = τs, w = w,
                samplesize = samplesize, K = K, metric = metric, α = α,
                p = p, undersampling = false)
    ε★s[:,counter] = ε★

    # zero-padding of ⟨ε★⟩ in order to also cover τ=0 (important for the multivariate case)
    ε★ = vec([0; ε★])
    # get the L-statistic-decrease for each peak in ⟨ε★⟩ and take the maximum one
    L_trials, max_idx = local_L_statistics(ε★, Y_act, s, τs, KNN, w, samplesize, metric)
    L_min, min_idx = findmin(L_trials)

    push!(τ_vals, τs[max_idx[min_idx]-1])
    push!(ts_vals, 1)
    push!(Ls, L_min)

    # create phase space vector for this embedding cycle
    Y_act = DelayEmbeddings.hcat_lagged_values(Y_act, s, τ_vals[counter+1])

    return Y_act
end

"""
Perform one embedding cycle on `Y_act` with a multivariate set Ys
"""
function pecuzal_multivariate_embedding_cycle!(
        Y_act, flag, Ys, τs, w, counter, ε★s, τ_vals, metric,
        Ls, ts_vals, samplesize, K, α, p, KNN)

    M = size(Ys,2)
    # in the 1st cycle we have to check all (size(Y,2)^2 combinations and pick
    # the tau according to maximum L-statistic decrease)
    if counter == 1
        Y_act = first_embedding_cycle_pecuzal!(Ys, M, τs, w, samplesize, K,
                                metric, α, p, KNN, τ_vals, ts_vals, Ls, ε★s)
    # in all other cycles we just have to check (size(Y,2)) combinations and pick
    # the tau according to maximum L-statistic decrease
    else
        Y_act = embedding_cycle_pecuzal!(Y_act, Ys, counter, M, τs, w, samplesize,
                            K, metric, α, p, KNN, τ_vals, ts_vals, Ls, ε★s)
    end
    return Y_act
end

"""
Perform the first embedding cycle of the multivariate embedding. Return the
actual reconstruction vector `Y_act`.
"""
function first_embedding_cycle_pecuzal!(Ys, M, τs, w, samplesize, K,
                        metric, α, p, KNN, τ_vals, ts_vals, Ls, ε★s)
    counter = 1
    L_min = zeros(M)
    L_act = zeros(M)
    L_min_idx = zeros(Int, M)
    ε★ = zeros(length(τs), M*M)
    idx = zeros(Int, M)
    ξ_min = zeros(M)
    for ts = 1:M
        ε★[:,1+(M*(ts-1)):M*ts], _ = pecora(Ys, (0,), (ts,); delays = τs,
                    w = w, samplesize = samplesize, K = K, metric = metric,
                    α = α, p = p, undersampling = false)
        L_min[ts], L_min_idx[ts], idx[ts]  = choose_right_embedding_params(
                                        ε★[:,1+(M*(ts-1)):M*ts], Ys[:,ts],
                                        Ys, τs, KNN, w, samplesize,
                                        metric)
    end
    L_mini, min_idx = findmin(L_min)
    # update τ_vals, ts_vals, Ls, ε★s
    push!(τ_vals, τs[L_min_idx[min_idx]])
    push!(ts_vals, min_idx)             # time series to start with
    push!(ts_vals, idx[min_idx])        # result of 1st embedding cycle
    push!(Ls, L_mini)                   # L-value of 1st embedding cycle
    ε★s[counter] = ε★[:,1+(M*(ts_vals[1]-1)):M*ts_vals[1]]

    # create phase space vector for this embedding cycle
    Y_act = DelayEmbeddings.hcat_lagged_values(Ys[:,ts_vals[counter]],
                                 Ys[:,ts_vals[counter+1]],τ_vals[counter+1])

    return Y_act
end

"""
Perform an embedding cycle of the multivariate embedding, but the first one.
Return the actual reconstruction vector `Y_act`.
"""
function embedding_cycle_pecuzal!(Y_act, Ys, counter, M, τs, w, samplesize,
                    K, metric, α, p, KNN, τ_vals, ts_vals, Ls, ε★s)

    ε★, _ = pecora(Ys, Tuple(τ_vals), Tuple(ts_vals); delays = τs, w = w,
            samplesize = samplesize, K = K, metric = metric, α = α,
            p = p, undersampling = false)
    # update τ_vals, ts_vals, Ls, ε★s
    choose_right_embedding_params!(ε★, Y_act, Ys, τ_vals, ts_vals, Ls, ε★s,
                                counter, τs, KNN, w, samplesize, metric)
    # create phase space vector for this embedding cycle
    Y_act = DelayEmbeddings.hcat_lagged_values(Y_act, Ys[:, ts_vals[counter+1]],
                                                        τ_vals[counter+1])
    return Y_act
end


"""
    Choose the maximum L-decrease and corresponding τ for each ε★-statistic,
based on picking the peak in ε★, which corresponds to the minimal `L`-statistic.
"""
function choose_right_embedding_params!(ε★, Y_act, Ys, τ_vals, ts_vals, Ls, ε★s,
                                 counter, τs, KNN, w, samplesize, metric)
    L_min_ = zeros(size(Ys,2))
    τ_idx = zeros(Int,size(Ys,2))
    for ts = 1:size(Ys,2)
        # zero-padding of ⟨ε★⟩ in order to also cover τ=0 (important for the multivariate case)
        # get the L-statistic for each peak in ⟨ε★⟩ and take the one according to L_min
        L_trials_, max_idx_ = local_L_statistics(vec([0; ε★[:,ts]]), Y_act, Ys[:,ts],
                                        τs, KNN, w, samplesize, metric)
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

"""
    Choose the right embedding parameters of the ε★-statistic in the first
embedding cycle. Return the `L`-decrease-value, the corresponding index value of
the chosen peak `τ_idx` and the number of the chosen time series to start with
`idx`.
"""
function choose_right_embedding_params(ε★, Y_act, Ys, τs, KNN, w, samplesize, metric)
    L_min_ = zeros(size(Ys,2))
    τ_idx = zeros(Int,size(Ys,2))
    for ts = 1:size(Ys,2)
        # zero-padding of ⟨ε★⟩ in order to also cover τ=0 (important for the multivariate case)
        # get the L-statistic for each peak in ⟨ε★⟩ and take the one according to L_min
        L_trials_, max_idx_ = local_L_statistics(vec([0; ε★[:,ts]]), Dataset(Y_act), Ys[:,ts],
                                        τs, KNN, w, samplesize, metric)
        L_min_[ts], min_idx_ = findmin(L_trials_)
        τ_idx[ts] = max_idx_[min_idx_]-1
    end
    idx = sortperm(L_min_)
    return L_min_[idx[1]], τ_idx[idx[1]], idx[1]
end


"""
    Return the maximum decrease of the L-statistic `L_decrease` and corresponding
delay-indices `max_idx` for all local maxima in ε★
"""
function local_L_statistics(ε★, Y_act, s, τs, KNN, w, samplesize, metric)
    maxima, max_idx = get_maxima(ε★) # determine local maxima in ⟨ε★⟩
    L_decrease = zeros(Float64, length(max_idx))
    for (i,τ_idx) in enumerate(max_idx)
        # create candidate phase space vector for this peak/τ-value
        Y_trial = DelayEmbeddings.hcat_lagged_values(Y_act, s, τs[τ_idx-1])
        # compute L-statistic for Y_act and Y_trial and get the maximum decrease
        L_decrease[i] = get_minimum_L_by_separation(Y_act, Y_trial, τs; K = KNN,
                                w = w, samplesize = samplesize, metric = metric)
    end
    return L_decrease, max_idx
end

"""
    Compute the `L`-Statistics for trajectories `Y_act` and `Y_trial` for many
prediction time horizons `Tw`. Return the (first) largest difference
L_trial - L_act.
"""
function get_minimum_L_by_separation(Y_act::Dataset{D, T}, Y_trial::Dataset{D2, T}, τs;
                    K::Int = 3, w::Int = 1, samplesize::Real = 1,
                    metric = Euclidean()) where {D, D2, T<:Real}
    # loop over time horizons until the maximum L-separation is reached
    L1, L2, dist, dist_former = 0, 0, 0, 99999999
    for Tw = 2:τs[end]
        L1 = uzal_cost(Y_act; Tw = Tw, K = K, w = w,
                samplesize = samplesize, metric = metric)
        L2 = uzal_cost(Y_trial; Tw = Tw, K = K, w = w,
                samplesize = samplesize, metric = metric)
        dist = L2 - L1
        if dist > dist_former && dist_former<0
            break
        else
            dist_former = dist
        end
    end
    return dist_former
end


function pecuzal_break_criterion(Ls, counter, max_num_of_cycles)
    flag = true
    if counter == 1 && Ls[end] > 0
        println("Algorithm stopped due to increasing L-values. "*
                "Valid embedding NOT achieved ⨉.")
        flag = false
    end
    if counter > 1 && Ls[end] > 0
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
Return the maxima of the given time series s and its indices
"""
function get_maxima(s::Vector{T}) where {T}
    maximas = T[]
    maximas_idx = Int[]
    N = length(s)
    flag = false
    first_point = 0
    for i = 2:N-1
        if s[i-1] < s[i] && s[i+1] < s[i]
            flag = false
            push!(maximas, s[i])
            push!(maximas_idx, i)
        end
        # handling constant values
        if flag
            if s[i+1] < s[first_point]
                flag = false
                push!(maximas, s[first_point])
                push!(maximas_idx, first_point)
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
        maximas, maximas_idx = findmax(s)
    end
    return maximas, maximas_idx
end


# """
#     uzal_cost(Y::Dataset; kwargs...) → L
# Compute the L-statistic `L` for input dataset `Y` according to Uzal et al.[^Uzal2011], based on
# theoretical arguments on noise amplification, the complexity of the
# reconstructed attractor and a direct measure of local stretch which constitutes
# an irrelevance measure. It serves as a cost function of a state space
# trajectory/embedding and therefore allows to estimate a "goodness of a
# embedding" and also to choose proper embedding parameters, while minimizing
# `L` over the parameter space. For receiving the local cost function `L_local`
# (for each point in state space - not averaged), use `uzal_cost_local(...)`.
#
# ## Keyword arguments
#
# * `samplesize = 0.5`: Number of considered fiducial points v as a fraction of
#   input state space trajectory `Y`'s length, in order to average the conditional
#   variances and neighborhood sizes (read algorithm description) to produce `L`.
# * `K = 3`: the amount of nearest neighbors considered, in order to compute σ_k^2
#   (read algorithm description).
#   If given a vector, minimum result over all `k ∈ K` is returned.
# * `metric = Euclidean()`: metric used for finding nearest neigbhors in the input
#   state space trajectory `Y.
# * `w = 1`: Theiler window (neighbors in time with index `w` close to the point,
#   that are excluded from being true neighbors). `w=0` means to exclude only the
#   point itself, and no temporal neighbors.
# * `Tw = 40`: The time horizon (in sampling units) up to which E_k^2 gets computed
#   and averaged over (read algorithm description).
#
# ## Description
# The `L`-statistic is based on theoretical arguments on noise amplification, the
# complexity of the reconstructed attractor and a direct measure of local stretch
# which constitutes an irrelevance measure. Technically, it is the logarithm of
# the product of `σ`-statistic and a normalization statistic `α`:
#
# L = log10(σ*α)
#
# The `σ`-statistic is computed as follows. `σ = √σ² = √(E²/ϵ²)`.
# `E²` approximates the conditional variance at each point in state space and
# for a time horizon `T ∈ Tw`, using `K` nearest neighbors. For each reference
# point of the state space trajectory, the neighborhood consists of the reference
# point itself and its `K+1` nearest neighbors. `E²` measures how strong
# a neighborhood expands during `T` time steps. `E²` is averaged over many time
# horizons `T = 1:Tw`. Consequently, `ϵ²` is the size of the neighborhood at the
# reference point itself and is defined as the mean pairwise distance of the
# neighborhood. Finally, `σ²` gets averaged over a range of reference points on
# the attractor, which is controlled by `samplesize`. This is just for performance
# reasons and the most accurate result will obviously be gained when setting
# `samplesize=1.0`
#
# The `α`-statistic is a normalization factor, such that `σ`'s from different
# embeddings can be compared. `α²` is defined as the inverse of the sum of
# the inverse of all `ϵ²`'s for all considered reference points.
#
# [^Uzal2011]: Uzal, L. C., Grinblat, G. L., Verdes, P. F. (2011). [Optimal reconstruction of dynamical systems: A noise amplification approach. Physical Review E 84, 016223](https://doi.org/10.1103/PhysRevE.84.016223).
# """
# function uzal_cost(Y::Dataset{D, ET};
#         Tw::Int = 40, K::Int = 3, w::Int = 1, samplesize::Real = 0.5,
#         metric = Euclidean()
#     ) where {D, ET}
#
#     # select a random state space vector sample according to input samplesize
#     NN = length(Y)-Tw;
#     NNN = floor(Int, samplesize*NN)
#     ns = sample(1:NN, NNN; replace=false) # the fiducial point indices
#
#     vs = Y[ns] # the fiducial points in the data set
#
#     vtree = KDTree(Y[1:end-Tw], metric)
#     allNNidxs, allNNdist = all_neighbors(vtree, vs, ns, K, w)
#     ϵ² = zeros(NNN)             # neighborhood size
#     E²_avrg = zeros(NNN)        # averaged conditional variance
#     E² = zeros(Tw)
#     ϵ_ball = zeros(ET, K+1, D) # preallocation
#     u_k = zeros(ET, D)
#
#     # loop over each fiducial point
#     for (i,v) in enumerate(vs)
#         NNidxs = allNNidxs[i] # indices of k nearest neighbors to v
#         # pairwise distance of fiducial points and `v`
#         pdsqrd = fiducial_pairwise_dist_sqrd(view(Y.data, NNidxs), v, metric)
#         ϵ²[i] = (2/(K*(K+1))) * pdsqrd  # Eq. 16
#         # loop over the different time horizons
#         for T = 1:Tw
#             E²[T] = comp_Ek2!(ϵ_ball, u_k, Y, ns[i], NNidxs, T, K, metric) # Eqs. 13 & 14
#         end
#         # Average E²[T] over all prediction horizons
#         E²_avrg[i] = mean(E²)                   # Eq. 15
#     end
#     σ² = E²_avrg ./ ϵ² # noise amplification σ², Eq. 17
#     σ²_avrg = mean(σ²) # averaged value of the noise amplification, Eq. 18
#     α² = 1 / sum(ϵ².^(-1)) # for normalization, Eq. 21
#     L = log10(sqrt(σ²_avrg)*sqrt(α²))
# end
#
# function fiducial_pairwise_dist_sqrd(fiducials, v, metric)
#     pd = zero(eltype(fiducials[1]))
#     pd += evaluate(metric, v, v)^2
#     for (i, v1) in enumerate(fiducials)
#         pd += evaluate(metric, v1, v)^2
#         for j in i+1:length(fiducials)
#             @inbounds pd += evaluate(metric, v1, fiducials[j])^2
#         end
#     end
#     return sum(pd)
# end
#
# """
#     comp_Ek2!(ϵ_ball,u_k,Y,v,NNidxs,T,K,metric) → E²(T)
# Returns the approximated conditional variance for a specific point in state space
# `ns` (index value) with its `K`-nearest neighbors, which indices are stored in
# `NNidxs`, for a time horizon `T`. This corresponds to Eqs. 13 & 14 in [^Uzal2011].
# The norm specified in input parameter `metric` is used for distance computations.
# """
# function comp_Ek2!(ϵ_ball, u_k, Y, ns::Int, NNidxs, T::Int, K::Int, metric)
#     # determine neighborhood `T` time steps ahead
#     ϵ_ball[1, :] .= Y[ns+T]
#     @inbounds for (i, j) in enumerate(NNidxs)
#         ϵ_ball[i+1, :] .= Y[j + T]
#     end
#
#     # compute center of mass
#     @inbounds for i in 1:size(Y)[2]; u_k[i] = sum(view(ϵ_ball, :, i))/(K+1); end # Eq. 14
#
#     E²_sum = 0
#     @inbounds for j = 1:K+1
#         E²_sum += (evaluate(metric, view(ϵ_ball, j, :), u_k))^2
#     end
#     E² = E²_sum / (K+1)         # Eq. 13
# end
