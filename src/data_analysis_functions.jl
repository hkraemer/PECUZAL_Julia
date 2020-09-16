using DrWatson
@quickactivate "new-embedding-methods"

using DynamicalSystems
using RecurrenceAnalysis
using DelayEmbeddings


"""
Perform the RecurrenceAnalysis of some reconstruction trajectories `Y₁`, `Y₂`,
`Y₃`, `Y₄`. Specifically, compute the fraction of recurrence rates from the
"original"/reference trajectory `Y_ref` with the one from the JRP of the
original `Y₁`, `Y₂`, `Y₃`, `Y₄` together with the reconstructed trajectory. Also
compute the recurrence time entropy of the recurrence plots of `Y₁`, `Y₂`, `Y₃`,
`Y₄` and `Y_ref`.

Keyword arguments:
*`ε = 0.05`: The used threshold for constructing the recurrence plots.
           The reconstruction method is fixed recurrence rate.
*`w_ref = 1`: Theiler window for the Dataset corresponding to `Y_ref`.
*`w₁ = 1`: Theiler window for the Dataset corresponding to `Y₁`.
*`w₂ = 1`: Theiler window for the Dataset corresponding to `Y₂`.
*`w₃ = 1`: Theiler window for the Dataset corresponding to `Y₃`.
*`w₄ = 1`: Theiler window for the Dataset corresponding to `Y₄`.
*`lmin = 2`: Minimum used line length for digaonal line based RQA measures.
"""
function perform_recurrence_analysis(Y_ref::Dataset, Y₁::Dataset,
                        Y₂::Dataset, Y₃::Dataset, Y₄::Dataset;
                        ε::Real = 0.05, w_ref::Int = 1, w₁::Int = 1,  w₂::Int = 1,
                        w₃::Int = 1, w₄::Int = 1, lmin::Int = 2)
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

    RQA_ref = rqa(R_ref; theiler = w_ref, lmin = lmin)
    RQA1 = rqa(R1; theiler = w₁, lmin = lmin)
    RQA2 = rqa(R2; theiler = w₂, lmin = lmin)
    RQA3 = rqa(R3; theiler = w₃, lmin = lmin)
    RQA4 = rqa(R4; theiler = w₄, lmin = lmin)

    return f1, f2, f3, f4, RQA_ref, RQA1, RQA2, RQA3, RQA4
end


"""
    standard_embedding_hegger(s::Vector; kwargs...) → `Y`, `τ`
Compute the reconstructed trajectory from a time series using the standard time
delay embedding. The delay `τ` is taken as the 1st minimum of the mutual
information [`estimate_dimension`](@ref) and the embedding dimension `m` is
estimated by using an FNN method from [^Hegger1999] [`fnn_uniform_hegger`](@ref)
with an optional keyword `fnn_thres = 0.05`, which defines at which fraction of
FNNs the search should break.
Return the reconstructed trajectory `Y` and the delay `τ`.

[^Hegger1999]: Hegger, Rainer and Kantz, Holger (1999). [Improved false nearest neighbor method to detect determinism in time series data. Physical Review E 60, 4970](https://doi.org/10.1103/PhysRevE.60.4970).
"""
function standard_embedding_hegger(s::Vector{T}; fnn_thres::Real = 0.05) where {T}
    τ = estimate_delay(s, "mi_min")
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

[^Hegger1999]: Hegger, Rainer and Kantz, Holger (1999). [Improved false nearest neighbor method to detect determinism in time series data. Physical Review E 60, 4970](https://doi.org/10.1103/PhysRevE.60.4970).
"""
function standard_embedding_cao(s::Vector{T}; cao_thres::Real = 0.05, m_max::Int = 10) where {T}
    τ = estimate_delay(s, "mi_min")
    rat = estimate_dimension(s, τ, 1:m_max, "afnn")
    for i = 1:m_max
        if abs(1-rat[i]) < cao_thres
            global m = i
            break
        end
    end
    if m > 1
        Y = embed(s, m, τ)
    else
        Y = s
    end
    return Y, τ
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
Compute elementwise product of two Recurrence Plots (==JRP)
"""
function elementwise_product(RP₁::RecurrenceMatrix, RP₂::RecurrenceMatrix)
    @assert size(RP₁) == size(RP₂)
    JRP = zeros(Int, size(RP₁))
    @inbounds for j = 1:size(RP₁,2)
        @inbounds for i = 1:size(RP₁,1)
            JRP[i,j] = RP₁[i,j] * RP₂[i,j]
        end
    end
    return JRP
end
