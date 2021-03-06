using DrWatson
@quickactivate "PECUZAL_Julia"

# Core functions for the usage in `/scripts`-folder

using DynamicalSystems
using RecurrenceAnalysis
using DelayEmbeddings
using Distances
using Statistics
using LinearAlgebra


"""
Perform the RecurrenceAnalysis of some reconstruction trajectories `Y₁`, `Y₂`,
`Y₃`, `Y₄`. Specifically, compute the fraction of recurrence rates from the
"original"/reference trajectory `Y_ref` with the one from the JRP of the
original `Y₁`, `Y₂`, `Y₃`, `Y₄` together with the reconstructed trajectory. Also
compute the recurrence time entropy of the recurrence plots of `Y₁`, `Y₂`, `Y₃`,
`Y₄` and `Y_ref`.

Keyword arguments:
*`ε = 0.05`: The used threshold for constructing the recurrence plots
    The reconstruction method is fixed recurrence rate.
*`w = 1`: Theiler window used for all Datasets
*`lmin = 2`: Minimum used line length for digaonal line based RQA measures.
*`kNN = 1`: The number of nearest neighbors used for obtaining the mutual
    nearest neighbors measure
"""
function perform_recurrence_analysis(Y_ref::Dataset, Y₁::Dataset,
                        Y₂::Dataset, Y₃::Dataset, Y₄::Dataset;
                        ε::Real = 0.05, w::Int = 1, lmin::Int = 2, kNN::Int = 1)

    N1 = length(Y₁)
    N2 = length(Y₂)
    N3 = length(Y₃)
    N4 = length(Y₄)
    N = minimum(hcat(N1, N2, N3, N4))

    R_ref = RecurrenceMatrix(Y_ref[1:N,:], ε; metric = "euclidean", fixedrate = true)
    R1 = RecurrenceMatrix(Y₁[1:N,:], ε; fixedrate = true)
    R2 = RecurrenceMatrix(Y₂[1:N,:], ε; fixedrate = true)
    R3 = RecurrenceMatrix(Y₃[1:N,:], ε; fixedrate = true)
    R4 = RecurrenceMatrix(Y₄[1:N,:], ε; fixedrate = true)

    f1 = jrp_rr_frac(R_ref, R1)
    f2 = jrp_rr_frac(R_ref, R2)
    f3 = jrp_rr_frac(R_ref, R3)
    f4 = jrp_rr_frac(R_ref, R4)

    mfnn1 = mfnn(Y_ref[1:N,:], Y₁[1:N,:]; w = w, kNN = kNN)
    mfnn2 = mfnn(Y_ref[1:N,:], Y₂[1:N,:]; w = w, kNN = kNN)
    mfnn3 = mfnn(Y_ref[1:N,:], Y₃[1:N,:]; w = w, kNN = kNN)
    mfnn4 = mfnn(Y_ref[1:N,:], Y₄[1:N,:]; w = w, kNN = kNN)

    RQA_ref = rqa(R_ref; theiler = w, lmin = lmin)
    RQA1 = rqa(R1; theiler = w, lmin = lmin)
    RQA2 = rqa(R2; theiler = w, lmin = lmin)
    RQA3 = rqa(R3; theiler = w, lmin = lmin)
    RQA4 = rqa(R4; theiler = w, lmin = lmin)

    return mfnn1, mfnn2, mfnn3, mfnn4, f1, f2, f3, f4, RQA_ref, RQA1, RQA2,
                                            RQA3, RQA4, R_ref, R1, R2, R3, R4
end


"""
    standard_embedding_hegger(s::Vector; kwargs...) → `Y`, `τ`
Compute the reconstructed trajectory from a time series using the standard time
delay embedding. The delay `τ` is taken as the 1st minimum of the mutual
information [`estimate_dimension`](@ref) and the embedding dimension `m` is
estimated by using an FNN method from [^Hegger1999] [`fnn_uniform_hegger`](@ref).
Return the reconstructed trajectory `Y` and the delay `τ`.

Keyword arguments:

*`fnn_thres = 0.05`: a threshold defining at which fraction of FNNs the search
    should break.
* The `method` can be one of the following:
* `"ac_zero"` : first delay at which the auto-correlation function becomes <0.
* `"ac_min"` : delay of first minimum of the auto-correlation function.
* `"mi_min"` : delay of first minimum of mutual information of `s` with itself
  (shifted for various `τs`). <- Default

[^Hegger1999]: Hegger, Rainer and Kantz, Holger (1999). [Improved false nearest neighbor method to detect determinism in time series data. Physical Review E 60, 4970](https://doi.org/10.1103/PhysRevE.60.4970).
"""
function standard_embedding_hegger(s::Vector{T}; method::String = "mi_min",
                                fnn_thres::Real = 0.05, τs = 1:200) where {T}
    @assert method=="ac_zero" || method=="mi_min" || method=="ac_min"
    "The absolute correlation function has elements that are = 0. "*
    "We can't fit an exponential to it. Please choose another method."

    τ = estimate_delay(s, method, τs)
    _, _, Y = fnn_uniform_hegger(s, τ; fnn_thres = fnn_thres)
    return Y, τ
end


"""
    standard_embedding_cao(s::Vector; kwargs...) → `Y`, `τ`
Compute the reconstructed trajectory from a time series using the standard time
delay embedding. The delay `τ` is taken as the 1st minimum of the mutual
information [`estimate_dimension`](@ref) and the embedding dimension `m` is
estimated by using an FNN method from Cao [`estimate_dimension`](@ref), with the
threshold parameter `cao_thres`.
Return the reconstructed trajectory `Y` and the delay `τ`.

Keyword arguments:
*`cao_thres = 0.05`: This threshold determines the tolerable deviation of the
    proposed statistic from the optimal value of 1, for breaking the algorithm.
*`m_max = 10`: The maximum embedding dimension, which is encountered by the
    algorithm.
* The `method` can be one of the following:
* `"ac_zero"` : first delay at which the auto-correlation function becomes <0.
* `"ac_min"` : delay of first minimum of the auto-correlation function.
* `"mi_min"` : delay of first minimum of mutual information of `s` with itself
  (shifted for various `τs`). <- Default

[^Hegger1999]: Hegger, Rainer and Kantz, Holger (1999). [Improved false nearest neighbor method to detect determinism in time series data. Physical Review E 60, 4970](https://doi.org/10.1103/PhysRevE.60.4970).
"""
function standard_embedding_cao(s::Vector{T}; cao_thres::Real = 0.05,
             τs = 1:200, method::String = "mi_min", m_max::Int = 10) where {T}
    @assert method=="ac_zero" || method=="mi_min" || method=="ac_min"
    "The absolute correlation function has elements that are = 0. "*
    "We can't fit an exponential to it. Please choose another method."

    τ = estimate_delay(s, method, τs)
    rat = estimate_dimension(s, τ, 1:m_max, "afnn")
    for i = 1:m_max
        if abs(1-rat[i]) < cao_thres
            global m = i
            break
        end
    end
    try
        if m > 1
            global Y = embed(s, m, τ)
        else
            global Y = s
        end
    catch
        global Y = s
    end
    return Y, τ, rat
end


"""
    fnn_uniform_hegger(s::Vector, τ::Int; kwargs...) →  `m`, `FNNs`, `Y`
Compute and return the optimal embedding dimension `m` for the time series `s`
and a uniform time delay `τ` after [^Hegger1999]. Return the optimal `m` and the
corresponding reconstruction vector `Y` according to that `m` and the input `τ`.
The optimal `m` is chosen, when the fraction of `FNNs` falls below the threshold
`fnn_thres` or when fraction of FNN's increases.

Keyword argument:
*`fnn_thres = 0.05`: Threshold, which defines the tolerable fraction of FNN's
    for which the algorithm breaks.
*`max_dimension = 10`: The maximum dimension which is encountered by the
    algorithm and after which it breaks, if the breaking criterion has not been
    met yet.
*`r = 2`: Obligatory threshold, which determines the maximum tolerable spreading
    of trajectories in the reconstruction space.
*`metric = Euclidean`: The norm used for distance computations.
*`w = 1` = The Theiler window, which excludes temporally correlated points from
    the nearest neighbor search.

[^Hegger1999]: Hegger, Rainer and Kantz, Holger (1999). [Improved false nearest neighbor method to detect determinism in time series data. Physical Review E 60, 4970](https://doi.org/10.1103/PhysRevE.60.4970).
"""
function fnn_uniform_hegger(s::Vector{T}, τ::Int; max_dimension::Int = 10,
            r::Real = 2, w::Int = 1, fnn_thres::Real = 0.05, metric = Euclidean()) where {T}
    @assert max_dimension > 0
    s = (s .- mean(s)) ./ std(s)
    Y_act = s

    vtree = KDTree(Dataset(s), metric)
    _, NNdist_old = DelayEmbeddings.all_neighbors(vtree, Dataset(s), 1:length(s), 1, w)

    FNNs = zeros(max_dimension)
    for m = 2:max_dimension+1
        Y_act = DelayEmbeddings.hcat_lagged_values(Y_act, s, m*τ)
        Y_act = regularize(Y_act)
        vtree = KDTree(Y_act, metric)
        _, NNdist_new = DelayEmbeddings.all_neighbors(vtree, Y_act, 1:length(Y_act), 1, w)

        FNNs[m-1] = DelayEmbeddings.fnn_embedding_cycle(view(NNdist_old,
                                            1:length(Y_act)), NNdist_new, r)

        flag = fnn_break_criterion(FNNs[1:m-1], fnn_thres)
        if flag
            global bm = m
            break
        else
            global bm = m
        end

        NNdist_old = NNdist_new
    end

    if bm>2
        Y_final = embed(s, bm-1, τ)
    else
        Y_final = s
    end
    return bm, FNNs[1:bm-1], Y_final
end

"""
Determines the break criterion for the Hegger-FNN-estimation
"""
function fnn_break_criterion(FNNs, fnn_thres)
    flag = false
    if FNNs[end] ≤ fnn_thres
        flag = true
        println("Algorithm stopped due to sufficiently small FNNs. "*
                "Valid embedding achieved ✓.")
    end
    if length(FNNs) > 1 && FNNs[end] > FNNs[end-1]
        flag = true
        println("Algorithm stopped due to rising FNNs. "*
                "Valid embedding achieved ✓.")
    end
    return flag
end


"""
Computes the similarity between recurrence plots `RP₁` and `RP₂`. Outputs the
fraction of recurrences rates gained from RP₁ and of the joint recurrence
plot `RP₁ .* RP₂`.
"""
function jrp_rr_frac(RP₁::RecurrenceMatrix, RP₂::RecurrenceMatrix)
    @assert size(RP₁) == size(RP₂)

    RR1 = sum(RP₁)/(size(RP₁,1)*size(RP₁,1))
    JRP = elementwise_product(RP₁, RP₂)
    RR2 = sum(JRP)/(size(JRP,1)*size(JRP,1))

    f = RR2 / RR1
    return f
end

"""
Computes the mututal false nearest neighbours (mfnn) for a reference trajectory
`Y_ref` and a reconstruction `Y_rec` after [^Rulkov1995].

Keyword arguments:

*`w = 1`: Theiler window for the surpression of serially correlated neighbors in
    the nearest neighbor-search
*`kNN = 1`: The number of considered nearest neighbours (in the paper always 1)

[^Rulkov1995]: Rulkov, Nikolai F. and Sushchik, Mikhail M. and Tsimring, Lev S. and Abarbanel, Henry D.I. (1995). [Generalized synchronization of chaos in directionally coupled chaotic systems. Physical Review E 51, 980](https://doi.org/10.1103/PhysRevE.51.980).
"""
function mfnn(Y_ref::Dataset, Y_rec::Dataset; w::Int = 1, kNN::Int = 1)

    @assert length(Y_ref) == length(Y_rec)
    @assert kNN > 0
    N = length(Y_ref)
    metric = Euclidean()

    # compute nearest neighbor distances for both trajectories
    vtree = KDTree(Y_ref, metric)
    allNNidxs_ref, _ = DelayEmbeddings.all_neighbors(vtree, Y_ref,
                                                        1:length(Y_ref), kNN, w)
    vtree = KDTree(Y_rec, metric)
    allNNidxs_rec, _ = DelayEmbeddings.all_neighbors(vtree, Y_rec,
                                                        1:length(Y_rec), kNN, w)

    F = zeros(N)
    factor1_nom = zeros(kNN)
    factor1_denom = zeros(kNN)
    factor2_nom = zeros(kNN)
    factor2_denom = zeros(kNN)
    for i = 1:N
        for j = 1:kNN
            factor1_nom[j] = evaluate(Euclidean(), Y_rec[i], Y_rec[allNNidxs_ref[i][j]])
            factor1_denom[j] = evaluate(Euclidean(), Y_ref[i], Y_ref[allNNidxs_ref[i][j]])
            factor2_nom[j] = evaluate(Euclidean(), Y_ref[i], Y_ref[allNNidxs_rec[i][j]])
            factor2_denom[j] = evaluate(Euclidean(), Y_rec[i], Y_rec[allNNidxs_rec[i][j]])
        end
        factor1 = sum(factor1_nom)/sum(factor1_denom)
        factor2 = sum(factor2_nom)/sum(factor2_denom)
        F[i] = factor1*factor2                                         # Eq.(27)
    end
    return mean(F)
end


"""
Compute elementwise product of two Recurrence Plots (==JRP)
"""
elementwise_product(RP₁, RP₂) = JointRecurrenceMatrix(RP₁, RP₂)

"""
Generate data from a AR(1) process for a initial value `u0`, a AR-coefficient
`α` and a white noise scaling parameter `p`. Return a time series of length `N`.
"""
function ar_process(u0::T, α::T, p::T, N::Int) where {T<:Real}
    x = zeros(T, N+10)
    x[1] = u0
    for i = 2:N+10
        x[i] = α*x[i-1] + p*randn()
    end
    return x[11:end]
end

"""
Generate data from a AR(2) process for a initial value `u0`, AR-coefficients
`α₁` and `α₂` and a white noise scaling parameter `p`. Return a time series of
length `N`.
"""
function ar_process2(u0::Vector, α₁::T, α₂::T, p::T, N::Int) where {T<:Real}
    @assert (length(u0) == 2)
    x = zeros(T, N+100)
    x[1] = u0[1]
    x[2] = u0[1]
    for i = 3:N+10
        x[i] = α₂*x[i-2] + α₁*x[i-1] + p*randn()
    end
    return x[101:end]
end

"""
Generate data from a AR(3) process for a initial value `u0`, AR-coefficients
`α₁`, `α₂` and `α₃` and a white noise scaling parameter `p`. Return a time series of
length `N`.
"""
function ar_process3(u0::Vector, α₁::T, α₂::T, α₃::T, p::T, N::Int) where {T<:Real}
    @assert (length(u0) == 3)
    x = zeros(T, N+100)
    x[1] = u0[1]
    x[2] = u0[2]
    x[3] = u0[3]
    for i = 4:N+10
        x[i] = α₃*x[i-3] + α₂*x[i-2] + α₁*x[i-1] + p*randn()
    end
    return x[101:end]
end

"""
Generate data from a AR(5) process for a initial value `u0`, AR-coefficients
`α₁`, `α₂`, `α₃`, `α₄`, `α₅` and a white noise scaling parameter
`p`. Return a time series of
length `N`.
"""
function ar_process5(u0::Vector, α₁::T, α₂::T, α₃::T, α₄::T, α₅::T, p::T, N::Int) where {T<:Real}
    @assert (length(u0) == 5)
    x = zeros(T, N+100)
    x[1] = u0[1]
    x[2] = u0[2]
    x[3] = u0[3]
    x[4] = u0[4]
    x[5] = u0[5]
    for i = 6:N+10
        x[i] = α₅*x[i-5] + α₄*x[i-4] + α₃*x[i-3] + α₂*x[i-2] + α₁*x[i-1] + p*randn()
    end
    return x[101:end]
end

"""
Generate data from a AR(7) process for a initial value `u0`, AR-coefficients
`α₁`, `α₂`, `α₃`, `α₄`, `α₅`, `α₆` and `α₇` and a white noise scaling parameter
`p`. Return a time series of
length `N`.
"""
function ar_process7(u0::Vector, α₁::T, α₂::T, α₃::T, α₄::T, α₅::T, α₆::T, α₇::T, p::T, N::Int) where {T<:Real}
    @assert (length(u0) == 7)
    x = zeros(T, N+100)
    x[1] = u0[1]
    x[2] = u0[2]
    x[3] = u0[3]
    x[4] = u0[4]
    x[5] = u0[5]
    x[6] = u0[6]
    x[6] = u0[7]
    for i = 8:N+10
        x[i] = α₇*x[i-7] + α₆*x[i-6] + α₅*x[i-5] + α₄*x[i-4] + α₃*x[i-3] + α₂*x[i-2] + α₁*x[i-1] + p*randn()
    end
    return x[101:end]
end

"""
compute_delta_L(s, τs, (js,) T_max; KNN = 3, w = 1, metric = Euclidean) → ΔL

    Compute the overall L-decrease `ΔL` of a given embedding of the time series
`s::Vector` with the delay values `τs` up to a maximum `T`-value `T_max`. We
respect the Theiler window `w`, the chosen `metric` and the number of considered
nearest neighbors `KNN`. It is also possible to compute `ΔL` for a multivariate
input `Y::Dataset`. Then one additionally needs to supply a vector `js`, which
lists the chosen time series corresponding to the given delay values in `τs`.
This is similar to the procedure in [`genembed`]@ref. The computations are based
on z-standardized input for ensuring comparability.
"""
function compute_delta_L(s::Vector{T}, τs::Vector{Int}, T_max::Int;
        KNN::Int = 3, w::Int = 1, metric = Euclidean()
    ) where {T}
    ts = regularize(Dataset(s))
    tss = vec(Matrix(Dataset(ts)))
    ΔL = 0
    Y_act = ts
    for i = 1 : (length(τs)-1)
        τ = τs[i+1]
        # embedding one cycle
        Y_next = DelayEmbeddings.hcat_lagged_values(Y_act, tss, τ)
        # compute ΔL for this cycle
        ΔL += uzal_cost_pecuzal(Y_act, Y_next, T_max; K = KNN, w = w, metric = metric)
        Y_act = Y_next
    end
    return ΔL
end

function compute_delta_L(Y::Dataset{D, T}, τs::Vector{Int}, js::Vector{Int},
        T_max::Int; KNN::Int = 3, w::Int = 1, metric = Euclidean()) where {D, T}
    @assert length(τs) == length(js)
    ts = regularize(Y)
    ΔL = 0
    Y_act = Dataset(ts[:,js[1]])
    for i = 1: (length(τs)-1)
        τ = τs[i+1]
        # embedding one cycle
        Y_next = DelayEmbeddings.hcat_lagged_values(Y_act, ts[:,js[i+1]], τ)
        # compute ΔL for this cycle
        ΔL += uzal_cost_pecuzal(Y_act, Y_next, T_max; K = KNN, w = w,
        metric = metric)
        Y_act = Y_next
    end
    return ΔL
end
