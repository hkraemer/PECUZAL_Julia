using DrWatson
@quickactivate "PECUZAL_Julia"

using Statistics
using DelimitedFiles

## compute mean and std of measures computed in `comm_roessler_ensemble_uni.jl`

# In the following the indexing for lines i and columns j is as follows:
# i = 1: σ = 0
# i = 2: σ = .1
# j = 1: TDE
# j = 2: Garcia & Almeida
# j = 3: MDOP
# j = 4: PECUZAL


# load variables
jrrf_tde = readdlm("./scripts/Roessler/ensemble results uni/jrrf_tde.csv")
jrrf_GA = readdlm("./scripts/Roessler/ensemble results uni/jrrf_GA.csv")
jrrf_mdop = readdlm("./scripts/Roessler/ensemble results uni/jrrf_mdop.csv")
jrrf_pec = readdlm("./scripts/Roessler/ensemble results uni/jrrf_pec.csv")
mfnn_tde = readdlm("./scripts/Roessler/ensemble results uni/mfnn_tde.csv")
mfnn_GA = readdlm("./scripts/Roessler/ensemble results uni/mfnn_GA.csv")
mfnn_mdop = readdlm("./scripts/Roessler/ensemble results uni/mfnn_mdop.csv")
mfnn_pec = readdlm("./scripts/Roessler/ensemble results uni/mfnn_pec.csv")
L_tde = readdlm("./scripts/Roessler/ensemble results uni/L_tde.csv")
L_GA = readdlm("./scripts/Roessler/ensemble results uni/L_GA.csv")
L_mdop = readdlm("./scripts/Roessler/ensemble results uni/L_mdop.csv")
L_pec = readdlm("./scripts/Roessler/ensemble results uni/L_pec.csv")
DET_tde = readdlm("./scripts/Roessler/ensemble results uni/DET_tde.csv")
DET_GA = readdlm("./scripts/Roessler/ensemble results uni/DET_GA.csv")
DET_mdop = readdlm("./scripts/Roessler/ensemble results uni/DET_mdop.csv")
DET_pec = readdlm("./scripts/Roessler/ensemble results uni/DET_pec.csv")
DET_ref = readdlm("./scripts/Roessler/ensemble results uni/DET_ref.csv")
LAM_tde = readdlm("./scripts/Roessler/ensemble results uni/LAM_tde.csv")
LAM_GA = readdlm("./scripts/Roessler/ensemble results uni/LAM_GA.csv")
LAM_mdop = readdlm("./scripts/Roessler/ensemble results uni/LAM_mdop.csv")
LAM_pec = readdlm("./scripts/Roessler/ensemble results uni/LAM_pec.csv")
LAM_ref = readdlm("./scripts/Roessler/ensemble results uni/LAM_ref.csv")
RTE_tde = readdlm("./scripts/Roessler/ensemble results uni/RTE_tde.csv")
RTE_GA = readdlm("./scripts/Roessler/ensemble results uni/RTE_GA.csv")
RTE_mdop = readdlm("./scripts/Roessler/ensemble results uni/RTE_mdop.csv")
RTE_pec = readdlm("./scripts/Roessler/ensemble results uni/RTE_pec.csv")
RTE_ref = readdlm("./scripts/Roessler/ensemble results uni/RTE_ref.csv")
ENTR_tde = readdlm("./scripts/Roessler/ensemble results uni/ENTR_tde.csv")
ENTR_GA = readdlm("./scripts/Roessler/ensemble results uni/ENTR_GA.csv")
ENTR_mdop = readdlm("./scripts/Roessler/ensemble results uni/ENTR_mdop.csv")
ENTR_pec = readdlm("./scripts/Roessler/ensemble results uni/ENTR_pec.csv")
ENTR_ref = readdlm("./scripts/Roessler/ensemble results uni/ENTR_ref.csv")
dims_tde = readdlm("./scripts/Roessler/ensemble results uni/dim_tde.csv")
dims_GA = readdlm("./scripts/Roessler/ensemble results uni/dim_GA.csv")
dims_mdop = readdlm("./scripts/Roessler/ensemble results uni/dim_mdop.csv")
dims_pec = readdlm("./scripts/Roessler/ensemble results uni/dim_pec.csv")
TRANS_tde = readdlm("./scripts/Roessler/ensemble results uni/TRANS_tde.csv")
TRANS_GA = readdlm("./scripts/Roessler/ensemble results uni/TRANS_GA.csv")
TRANS_mdop = readdlm("./scripts/Roessler/ensemble results uni/TRANS_mdop.csv")
TRANS_pec = readdlm("./scripts/Roessler/ensemble results uni/TRANS_pec.csv")
TRANS_ref = readdlm("./scripts/Roessler/ensemble results uni/TRANS_ref.csv")


mean_JRRF = zeros(2,4)
std_JRRF = zeros(2,4)
mean_mfnn = zeros(2,4)
std_mfnn = zeros(2,4)
mean_L = zeros(2,4)
std_L = zeros(2,4)
mean_dims = zeros(2,4)
std_dims = zeros(2,4)

for i = 1:4
	for j = 1:2
		if i == 1
			mean_JRRF[j,i] = mean(jrrf_tde[j,:])
			std_JRRF[j,i] = std(jrrf_tde[j,:])
			mean_mfnn[j,i] = mean(mfnn_tde[j,:])
			std_mfnn[j,i] = std(mfnn_tde[j,:])
			mean_L[j,i] = mean(L_tde[j,:])
			std_L[j,i] = std(L_tde[j,:])
			mean_dims[j,i] = mean(dims_tde[j,:])
			std_dims[j,i] = std(dims_tde[j,:])
		elseif i == 2
			mean_JRRF[j,i] = mean(jrrf_GA[j,:])
			std_JRRF[j,i] = std(jrrf_GA[j,:])
			mean_mfnn[j,i] = mean(mfnn_GA[j,:])
			std_mfnn[j,i] = std(mfnn_GA[j,:])
			mean_L[j,i] = mean(L_GA[j,:])
			std_L[j,i] = std(L_GA[j,:])
			mean_dims[j,i] = mean(dims_GA[j,:])
			std_dims[j,i] = std(dims_GA[j,:])
		elseif i == 3
			mean_JRRF[j,i] = mean(jrrf_mdop[j,:])
			std_JRRF[j,i] = std(jrrf_mdop[j,:])
			mean_mfnn[j,i] = mean(mfnn_mdop[j,:])
			std_mfnn[j,i] = std(mfnn_mdop[j,:])
			mean_L[j,i] = mean(L_mdop[j,:])
			std_L[j,i] = std(L_mdop[j,:])
			mean_dims[j,i] = mean(dims_mdop[j,:])
			std_dims[j,i] = std(dims_mdop[j,:])
		elseif i == 4
			mean_JRRF[j,i] = mean(jrrf_pec[j,:])
			std_JRRF[j,i] = std(jrrf_pec[j,:])
			mean_mfnn[j,i] = mean(mfnn_pec[j,:])
			std_mfnn[j,i] = std(mfnn_pec[j,:])
			mean_L[j,i] = mean(L_pec[j,:])
			std_L[j,i] = std(L_pec[j,:])
			mean_dims[j,i] = mean(dims_pec[j,:])
			std_dims[j,i] = std(dims_pec[j,:])
		end
	end
end

mean_DET = zeros(2,4)
std_DET = zeros(2,4)
mean_LAM = zeros(2,4)
std_LAM = zeros(2,4)
mean_RTE = zeros(2,4)
std_RTE = zeros(2,4)
mean_ENTR = zeros(2,4)
std_ENTR = zeros(2,4)
mean_TRANS = zeros(2,4)
std_TRANS  = zeros(2,4)

for i = 1:4
	for j = 1:2
		if i == 1
			mean_DET[j,i] = mean( abs.(DET_tde[j,:] .- DET_ref[j,:]) ./ DET_ref[j,:])
			std_DET[j,i] = std( abs.(DET_tde[j,:] .- DET_ref[j,:]) ./ DET_ref[j,:])
			mean_LAM[j,i] = mean( abs.(LAM_tde[j,:] .- LAM_ref[j,:]) ./ LAM_ref[j,:])
			std_LAM[j,i] = std( abs.(LAM_tde[j,:] .- LAM_ref[j,:]) ./ LAM_ref[j,:])
			mean_ENTR[j,i] = mean( abs.(ENTR_tde[j,:] .- ENTR_ref[j,:]) ./ ENTR_ref[j,:])
			std_ENTR[j,i] = std( abs.(ENTR_tde[j,:] .- ENTR_ref[j,:]) ./ ENTR_ref[j,:])
			mean_RTE[j,i] = mean( abs.(RTE_tde[j,:] .- RTE_ref[j,:]) ./ RTE_ref[j,:])
			std_RTE[j,i] = std( abs.(RTE_tde[j,:] .- RTE_ref[j,:]) ./ RTE_ref[j,:])
			mean_TRANS[j,i] = mean( abs.(TRANS_tde[j,:] .- TRANS_ref[j,:]) ./ TRANS_ref[j,:])
			std_TRANS[j,i] = std( abs.(TRANS_tde[j,:] .- TRANS_ref[j,:]) ./ TRANS_ref[j,:])
		elseif i == 2
			mean_DET[j,i] = mean( abs.(DET_GA[j,:] .- DET_ref[j,:]) ./ DET_ref[j,:])
			std_DET[j,i] = std( abs.(DET_GA[j,:] .- DET_ref[j,:]) ./ DET_ref[j,:])
			mean_LAM[j,i] = mean( abs.(LAM_GA[j,:] .- LAM_ref[j,:]) ./ LAM_ref[j,:])
			std_LAM[j,i] = std( abs.(LAM_GA[j,:] .- LAM_ref[j,:]) ./ LAM_ref[j,:])
			mean_ENTR[j,i] = mean( abs.(ENTR_GA[j,:] .- ENTR_ref[j,:]) ./ ENTR_ref[j,:])
			std_ENTR[j,i] = std( abs.(ENTR_GA[j,:] .- ENTR_ref[j,:]) ./ ENTR_ref[j,:])
			mean_RTE[j,i] = mean( abs.(RTE_GA[j,:] .- RTE_ref[j,:]) ./ RTE_ref[j,:])
			std_RTE[j,i] = std( abs.(RTE_GA[j,:] .- RTE_ref[j,:]) ./ RTE_ref[j,:])
			mean_TRANS[j,i] = mean( abs.(TRANS_GA[j,:] .- TRANS_ref[j,:]) ./ TRANS_ref[j,:])
			std_TRANS[j,i] = std( abs.(TRANS_GA[j,:] .- TRANS_ref[j,:]) ./ TRANS_ref[j,:])
		elseif i == 3
			mean_DET[j,i] = mean( abs.(DET_mdop[j,:] .- DET_ref[j,:]) ./ DET_ref[j,:])
			std_DET[j,i] = std( abs.(DET_mdop[j,:] .- DET_ref[j,:]) ./ DET_ref[j,:])
			mean_LAM[j,i] = mean( abs.(LAM_mdop[j,:] .- LAM_ref[j,:]) ./ LAM_ref[j,:])
			std_LAM[j,i] = std( abs.(LAM_mdop[j,:] .- LAM_ref[j,:]) ./ LAM_ref[j,:])
			mean_ENTR[j,i] = mean( abs.(ENTR_mdop[j,:] .- ENTR_ref[j,:]) ./ ENTR_ref[j,:])
			std_ENTR[j,i] = std( abs.(ENTR_mdop[j,:] .- ENTR_ref[j,:]) ./ ENTR_ref[j,:])
			mean_RTE[j,i] = mean( abs.(RTE_mdop[j,:] .- RTE_ref[j,:]) ./ RTE_ref[j,:])
			std_RTE[j,i] = std( abs.(RTE_mdop[j,:] .- RTE_ref[j,:]) ./ RTE_ref[j,:])
			mean_TRANS[j,i] = mean( abs.(TRANS_mdop[j,:] .- TRANS_ref[j,:]) ./ TRANS_ref[j,:])
			std_TRANS[j,i] = std( abs.(TRANS_mdop[j,:] .- TRANS_ref[j,:]) ./ TRANS_ref[j,:])
		elseif i == 4
			mean_DET[j,i] = mean( abs.(DET_pec[j,:] .- DET_ref[j,:]) ./ DET_ref[j,:])
			std_DET[j,i] = std( abs.(DET_pec[j,:] .- DET_ref[j,:]) ./ DET_ref[j,:])
			mean_LAM[j,i] = mean( abs.(LAM_pec[j,:] .- LAM_ref[j,:]) ./ LAM_ref[j,:])
			std_LAM[j,i] = std( abs.(LAM_pec[j,:] .- LAM_ref[j,:]) ./ LAM_ref[j,:])
			mean_ENTR[j,i] = mean( abs.(ENTR_pec[j,:] .- ENTR_ref[j,:]) ./ ENTR_ref[j,:])
			std_ENTR[j,i] = std( abs.(ENTR_pec[j,:] .- ENTR_ref[j,:]) ./ ENTR_ref[j,:])
			mean_RTE[j,i] = mean( abs.(RTE_pec[j,:] .- RTE_ref[j,:]) ./ RTE_ref[j,:])
			std_RTE[j,i] = std( abs.(RTE_pec[j,:] .- RTE_ref[j,:]) ./ RTE_ref[j,:])
			mean_TRANS[j,i] = mean( abs.(TRANS_pec[j,:] .- TRANS_ref[j,:]) ./ TRANS_ref[j,:])
			std_TRANS[j,i] = std( abs.(TRANS_pec[j,:] .- TRANS_ref[j,:]) ./ TRANS_ref[j,:])
		end
	end
end
