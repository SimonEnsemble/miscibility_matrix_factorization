### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ c3f6fd50-1697-11ed-1af7-c717d141dfac
using CSV, DataFrames,Statistics, LogExpFunctions, CairoMakie,LinearAlgebra, Random , StatsBase, PlutoUI, ColorSchemes, ScikitLearn

# ╔═╡ b683bdc4-9956-4176-9c03-bebf15d0243e
using Distributions

# ╔═╡ 57a6f607-cd48-4e57-a452-06e0fcdd435c
using KernelFunctions

# ╔═╡ 57504276-6927-496e-9441-79ed5fbdd78e
import AlgebraOfGraphics

# ╔═╡ 2cbce39d-e9c3-474a-9a0f-ceb6442c27de
AlgebraOfGraphics.set_aog_theme!(fonts=[AlgebraOfGraphics.firasans("Light"), AlgebraOfGraphics.firasans("Light")]); update_theme!(fontsize=16)

# ╔═╡ 7b533d68-4bb8-4b1d-a48f-0adffc99ae97
import MLJBase: partition, StratifiedCV, train_test_pairs

# ╔═╡ 9cf628c4-3194-423a-90fc-94bb46c25450
@sk_import metrics: (confusion_matrix, accuracy_score, f1_score, precision_score, recall_score)

# ╔═╡ 5df72253-8573-46e9-a21b-923e0ce3f360
@sk_import decomposition: PCA

# ╔═╡ cabd0924-cae4-4992-bf49-92353aeea388
TableOfContents()

# ╔═╡ adb50317-7042-4894-8f2a-e9f6505adcfe
label_to_color = Dict(
	0 => reverse(ColorSchemes.turbid)[1],
	1 => reverse(ColorSchemes.turbid)[end]
)

# ╔═╡ f924129b-4e85-45b4-95ce-6846d69d16f3
md"## read in, viz miscibility data
* `0` $\implies$ immisible
* `1` $\implies$ miscible.

[link to original paper](https://pubs.acs.org/doi/10.1021/acsami.0c21036).

`pairs.csv` is \"Supplemental Table 4 006: Sample Code (XLSX)\"

`compounds.csv` is from Supp Table 1 link (but improperly labeled as Supp Table 2...)
"

# ╔═╡ b9e723cf-3af9-40b7-a9de-61893c633d53
begin
	df = CSV.read("pairs.csv", DataFrame, header=0)
	M_complete = Int.(Array(df)) # miscibility matrix
end

# ╔═╡ e559a0c8-7e1f-4fd6-a029-35b6354cd609
n_compounds = size(M_complete)[1]

# ╔═╡ 61377868-d14c-4927-a6e1-6d887903fbf8
sum(M_complete .== 0)

# ╔═╡ 7b9c60c4-b7a7-499b-8888-508d3610a966
sum(M_complete .== 1)

# ╔═╡ 35c79c1f-f086-419f-a2dc-2c1710e3ba59
md"some tests..."

# ╔═╡ dc065d67-229c-4f3f-8912-45a2deb3bf5e
begin
	@assert M_complete' == M_complete   # symmetric
	@assert all(diag(M_complete) .== 1) # diagonal all ones
	# miscible entries
	@assert M_complete[2, 1] == M_complete[4, 1] == M_complete[13, 3] == 1
	# immiscible entries
	@assert M_complete[4, 3] == M_complete[13, 4] == M_complete[12, 3] == M_complete[15, 1] == 0
	@assert all(diag(M_complete) .== 1) # all compounds miscible with themselves.
end

# ╔═╡ 22531e96-f0f7-45c0-a327-53094a83ddcf
begin
	compounds = CSV.read("compounds.csv", DataFrame)

	# add ATPS (Aqueous Two-Phase System) forming fractions data
	compounds[:, :ATPS_frac] = 1 .- sum(M_complete, dims=2)[:] / n_compounds
	compounds
end

# ╔═╡ 01890165-4cb1-4917-959d-9a0c28e60073
begin
	feature_names = ["monomer_mw", "XlogP3", "h-bond_donors", 
					 "h-bond_acceptors", "complexity", "concentration",
					 "polymer_mw"]
	
	feature_matrix = Matrix(compounds[:, feature_names])
	# normalize
	for j = 1:size(feature_matrix)[2]
		μ = mean(feature_matrix[:, j])
		σ = std(feature_matrix[:, j])
		feature_matrix[:, j] .= (feature_matrix[:, j] .- μ) / σ
	end
	feature_matrix
end

# ╔═╡ 9b0388ac-a640-4ffa-8531-b1ae7b68393f
function viz_miscibility_matrix(M)
	hm_colormap = reverse(ColorSchemes.turbid)
	
	fig = Figure(resolution=(1200, 1200))
	ax  = Axis(fig[1, 1], 
		xlabel="compound", ylabel="compound", 
		title="miscibility matrix",
		xgridvisible=false, 
		ygridvisible=false, 
		aspect=DataAspect(),
		xticks=(1:n_compounds, compounds[:, "NAME"]),
		yticks=(1:n_compounds, reverse(compounds[:, "NAME"])),
		xticklabelrotation=π/2
	)
	heatmap!(reverse(M', dims=2), colormap=hm_colormap, aspect=DataAspect())
	for i = 1:n_compounds
		hlines!(i - 0.5, color="gray", linewidth=1)
		vlines!(i - 0.5, color="gray", linewidth=1)
	end
	## legend
	legend_patches = [
		PolyElement(color=hm_colormap[1], strokecolor="gray"), PolyElement(color=hm_colormap[end], strokecolor="gray")
	]
	legend_labels = ["immiscible", "miscible"]
	legend_patchsize = (35, 35)
	if any(ismissing.(M))
		push!(legend_patches, PolyElement(color="white", strokecolor="gray"))
		push!(legend_labels, "missing")
		legend_patchsize = (35, 35, 35)
	end
	Legend(fig[1, 2], legend_patches, legend_labels, patchsize=legend_patchsize)
	return fig
end

# ╔═╡ acf912f1-51cf-4943-8964-c51cc3bf3427
viz_miscibility_matrix(M_complete)

# ╔═╡ 3a6f1f39-b598-4ce2-8858-b17551133bfe
compound_names = String.(compounds[:, "NAME"])

# ╔═╡ 3153dd82-8be6-4683-b93d-9629c3d3312f
md"## similarity matrix"

# ╔═╡ d6d83535-fa87-4de9-8ddc-0920a6bad075
compound_classes = String.(compounds[:, "CLASS"])

# ╔═╡ e2e9ae6a-f156-4a91-a63c-6ab7365fe19b
unique(compound_classes)

# ╔═╡ e6608357-deeb-4bc4-b5b6-5b5b50f07a3d
function kernel_matrix(σ::Nothing)
	K = zeros(n_compounds, n_compounds)
	for i = 1:n_compounds
		for j = 1:n_compounds
			if i == j
				K[i, j] = NaN # for safety
			else
				K[i, j] = compound_classes[i] == compound_classes[j]
			end
		end
	end
	return K
end

# ╔═╡ 32837121-174a-498a-9586-a5e4ac13398c
function kernel_matrix(σ::Float64)
	kernel = σ*SqExponentialKernel() 
	
	
	K = zeros(n_compounds, n_compounds)
	for i = 1:n_compounds
		for j = 1:n_compounds
			if i == j
				K[i, j] = NaN # for safety
			else
				K[i, j] = kernel(feature_matrix[i, :]', feature_matrix[j, :]')
			end
		end
	end
	return K
end

# ╔═╡ 4aea5499-5cdc-41e4-ad84-ae0eaf6c69e5
kernel_matrix(nothing)

# ╔═╡ 1be7e102-d7e9-4c08-989b-211ae302754e
kernel_matrix(0.1)

# ╔═╡ ee675570-94ef-4061-9438-7e60a0e5dae1
md"## introducing missing values"

# ╔═╡ d1062c18-909e-49f7-b7b8-82a566316b64
mutable struct MiscibilityData
	M::Matrix{Union{Int64, Missing}}
	n_compounds::Int
	compound_names::Vector{String}
	compound_classes::Vector{String}
	X::Matrix{Float64} # feature matrix
	ids_obs::Vector{Tuple{Int, Int}}
	ids_missing::Vector{Tuple{Int, Int}}
	class_wt::Tuple{Float64, Float64}
end

# ╔═╡ 6bc8a921-999d-4be4-a35c-aa3a0383c53e
function check_consistency(data::MiscibilityData)
	@assert all(ismissing.([data.M[i, j] for (i, j) in data.ids_missing]))
	@assert all(ismissing.([data.M[j, i] for (i, j) in data.ids_missing]))
	
	@assert all(.! ismissing.([data.M[i, j] for (i, j) in data.ids_obs]))
	@assert all(.! ismissing.([data.M[j, i] for (i, j) in data.ids_obs]))
end

# ╔═╡ 885c20b3-6376-4ce3-991a-8e1db6e88771
# θ:fraction of missing entries
# important note: we only pay attn to upper triangle of matrix.
function sim_data_collection(θ::Float64, M_complete::Matrix{Int64}; seed::Int=97330)
	Random.seed!(seed)
	
	@assert (θ < 1.0) && (θ > 0.0)
	n_compounds = size(M_complete)[1]

	# make list of entries in upper diagonal of the matrix
	ids_all_pairs = [(i, j) for i = 1:n_compounds for j = i+1:n_compounds]
	@assert length(ids_all_pairs) == (n_compounds ^ 2 - n_compounds) ÷ 2

	# number of observed entries
	nb_observed = floor(Int, (1-θ) * length(ids_all_pairs))

	# sample observed tuples. stratified split.
	ms_all_pairs = [M_complete[i, j] for (i, j) in ids_all_pairs] # for stratifying
    ids_missing, ids_obs = partition(ids_all_pairs, θ, shuffle=true, 
		                             stratify=ms_all_pairs)
	
	# construct the matrix with missing values
    M = zeros(Union{Int64, Missing}, n_compounds, n_compounds)
	fill!(M, missing) # default missing
	# fill in observed values
	for (i, j) in ids_obs
		M[i, j] = M_complete[i, j]
		M[j, i] = M_complete[i, j]
	end
	# refill diagonal as miscible
	for i = 1:n_compounds
		M[i, i] = 1
	end

	# to down-weigh the majority class
	fraction_miscible = mean([M[i, j] for (i, j) in ids_obs])
	class_wt = (1.0, (1 - fraction_miscible) / fraction_miscible)
	
    return MiscibilityData(M, n_compounds, compound_names, compound_classes, feature_matrix, ids_obs, ids_missing, class_wt)
end

# ╔═╡ 35e3ff6c-7e08-4956-b6d2-e9c6b8b076c6
data = sim_data_collection(0.5, M_complete)

# ╔═╡ 8deb8f6f-2be0-4c48-bf2b-4bf2e391c5c1
any(ismissing.(data.M))

# ╔═╡ 4c217dfa-2f9c-4ac9-8dfa-3f6bf270139f
viz_miscibility_matrix(data.M)

# ╔═╡ 54b331da-d6c7-4b5d-b71a-bc7d33c3c71a
md"some tests"

# ╔═╡ c106b08f-2025-4531-90e2-ed1d4d8e8b36
begin
	@assert all(skipmissing(data.M') .== skipmissing(data.M))
	# @assert sum(ismissing.(data.M)) == Int((n_compounds + (n_compounds ^ 2 - n_compounds) / 2))
	for i = 1:n_compounds, j = 1:n_compounds
		if ! ismissing(data.M[i, j])
			@assert data.M[i, j] == M_complete[i, j]
		end
	end
	check_consistency(data)
end

# ╔═╡ eaa4560e-cb0a-47bc-87ba-9ff4e3faa2c2
md"## the loss function and its gradient"

# ╔═╡ fabe9de6-f0e3-49fd-84b5-fe973e39f5a0
mutable struct Model
	C::Matrix{Float64}
	b::Float64
	γ::Float64 # regularization param for graph structure
	λ::Float64 # regularization param on latent reps
	σ::Union{Float64, Nothing} # scale param
end

# ╔═╡ af50a2c4-d8aa-49bd-b596-39964d2d8c7e
function pred_mᵢⱼ(model::Model, i::Int, j::Int)
	cᵢ = model.C[:, i]
	cⱼ = model.C[:, j]
	return logistic(dot(cᵢ, cⱼ) + model.b)
end

# ╔═╡ 500a3554-0bbf-42d2-ae4b-0b6c59dbeea6
function loss(data::MiscibilityData, model::Model, K::Matrix{Float64})
	# code below assumes M all 0's or 1's
	@assert all(skipmissing(data.M) .== 1 .|| skipmissing(data.M) .== 0)
	# latent vectors are columns of C
	@assert size(model.C)[2] == size(data.M)[1] == size(data.M)[2]
	
	# cross-entropy loss expressing mismatch over observations
	l = 0.0
	for i = 1:data.n_compounds
		l += model.λ * norm(model.C[:, i])
		for j = (i+1):data.n_compounds
			###
			#    graph-regularization term
			#      (looks at ALL pairs)
			###
			cᵢ = model.C[:, i]
			cⱼ = model.C[:, j]
			l += model.γ * K[i, j] * sum((cᵢ - cⱼ) .^ 2)
			
			###
			#    mis-match term
			#      (looks at pair only if not missing.)
			###
			if ! ismissing(data.M[i, j])
				# model prediction
				m̂ᵢⱼ = pred_mᵢⱼ(model, i, j)
				
				# truth
				mᵢⱼ = data.M[i, j]
				
				# increment mismatch loss. wt by class.
				wt = data.class_wt[mᵢⱼ + 1]
				l -= (mᵢⱼ * log(m̂ᵢⱼ) + (1 - mᵢⱼ) * log(1 - m̂ᵢⱼ)) * wt
			end
		end
	end
	return l
end

# ╔═╡ 4758b6e9-3c01-42b4-85bf-30ac2298e511
function ∇_cᵢ(data::MiscibilityData, model::Model, c::Int, K::Matrix{Float64})
	∇ = zeros(size(model.C)[1])
	
	∇ += model.λ * 2 * model.C[:, c] # regularize latent rep.
	# cross-entropy loss expressing mismatch over observations
	for j = 1:data.n_compounds
		###
		#   contribution to gradient by reg. term 
		###
		if c != j
			∇ += model.γ * K[c, j] * 2 * (model.C[:, c] - model.C[:, j]) # graph
		end
		
		###
		#   contribution to gradient by mis-matches
		###
		if ! ismissing(data.M[c, j])
			# model prediction
			m̂ᵢⱼ = pred_mᵢⱼ(model, c, j)
			# truth
			mᵢⱼ = data.M[c, j]
			# increment loss
			cⱼ = model.C[:, j]
			∇ += cⱼ * (m̂ᵢⱼ - mᵢⱼ) * data.class_wt[mᵢⱼ + 1]
		end
	end
	
	return ∇
end

# ╔═╡ 81f11413-bc17-4d38-880b-1afbca4ee01d
function ∇_b(data::MiscibilityData, model::Model)
	∇ = 0.0
	
	# cross-entropy loss expressing mismatch over observations
	for i = 1:data.n_compounds
		for j = (i+1):data.n_compounds
			if ! ismissing(data.M[i, j])
				# model prediction
				m̂ᵢⱼ = pred_mᵢⱼ(model, i, j)
				# truth
				mᵢⱼ = data.M[i, j]
				# increment loss
				∇ += (m̂ᵢⱼ - mᵢⱼ) * data.class_wt[mᵢⱼ + 1]
			end
		end
	end
	return ∇
end

# ╔═╡ c44a5395-8143-4e33-9eff-dd7f01a1217b
function grad_descent_epoch!(data::MiscibilityData, model::Model, K::Matrix{Float64}; α::Float64=0.01)
	# update compound vectors
	for c = shuffle(1:data.n_compounds)
		model.C[:, c] -= α * ∇_cᵢ(data, model, c, K)
	end
	# update bias
	model.b -= α * ∇_b(data, model)
	return nothing
end

# ╔═╡ a87cc2ef-b91f-49a2-9fbe-cef50192976d
md"## training"

# ╔═╡ 8b9046d8-8a0b-4f4a-b461-8c20349ecd30
function fraction_miscible(M::Union{Matrix{Int}, Matrix{Union{Missing, Int}}})
	return sum(skipmissing(M) .== 1) / sum(.! ismissing.(M))
end

# ╔═╡ 52deabd5-bc18-4659-a7e5-9f42b1e5f591
function construct_train_model(hyperparams::NamedTuple, 
	                           data::MiscibilityData, 
	                           nb_epochs::Int; 
							   α::Float64=0.01)
	# initialize model
	f_miscible = fraction_miscible(data.M)
	b_guess = log(f_miscible / (1 - f_miscible))
	C_guess = 0.5 * rand(Uniform(-1, 1), (hyperparams.k, data.n_compounds))
	model = Model(C_guess, b_guess, hyperparams.γ, hyperparams.λ, hyperparams.σ)

	# compute kernel matrix. note σ can be nothing.
	K = kernel_matrix(hyperparams.σ)

	# gradient descent epochs. keep track of loss.
	losses = [NaN for _ = 1:nb_epochs]
	for s = 1:nb_epochs
		grad_descent_epoch!(data, model, K, α=α)
		losses[s] = loss(data, model, K)
	end
	return model, losses
end

# ╔═╡ 8db228a7-bf2b-4f1e-ab7e-596b2957498f
hyperparams = (k=2, γ=0.01, λ=.5, σ=0.1)

# ╔═╡ 660f9ba7-871a-43df-b69d-d792a8e75456
nb_epochs = 1000

# ╔═╡ ba7a3de9-b37f-4c67-869b-fceb7915ffab
model, losses = construct_train_model(hyperparams, data, nb_epochs) # just an example

# ╔═╡ 1a178491-815e-49e3-b423-832a24abdba9
function _viz_loss!(ax, losses::Vector{Float64})
	lines!(ax, 1:length(losses), losses)
end

# ╔═╡ 3f1084c4-28f8-4f0d-98b8-2a0426c89e18
function viz_loss(losses::Vector{Float64})
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="# epochs", ylabel="loss")
	_viz_loss!(ax, losses)
	fig
end

# ╔═╡ 0ede5258-7b6a-4975-8d6b-65639bb7ac74
viz_loss(losses)

# ╔═╡ c939cb70-858d-4b57-977a-693097c78499
md"## hyperparam optimization"

# ╔═╡ 6ceaf279-5e77-4585-812b-506c152fd6ab
function compute_perf_metrics(model::Model, 
	                          data::MiscibilityData, 
	                          ids::Vector{Tuple{Int, Int}})
	m = [M_complete[i, j]            for (i, j) in ids]
	m̂ = [pred_mᵢⱼ(model, i, j) > 0.5 for (i, j) in ids]

	return (f1=f1_score(m, m̂),
		    accuracy=accuracy_score(m, m̂),
	        precision=precision_score(m, m̂),
		    recall=recall_score(m, m̂)
	)
end

# ╔═╡ 53349731-ff98-4ed7-9877-fe8cf389e6e3
compute_perf_metrics(model, data, data.ids_missing)[:f1]

# ╔═╡ 5efedaad-892f-4005-8ea2-27ce5b3fd7d8
perf_metrics

# ╔═╡ 9eef7d3f-a121-41ff-828b-042910b56789
ngrid = 10

# ╔═╡ 92cd854b-fb84-4981-ba8e-b24ecd1509c7
function do_hyperparam_optimization(
	data::MiscibilityData,
	cv_hyperparams::Vector{<:NamedTuple};
	nfolds::Int=3,
	type_of_perf_metric::Symbol=:f1
)
	# cross-validation split
	kf_split = train_test_pairs(StratifiedCV(nfolds=nfolds, shuffle=true), 
		                    1:length(data.ids_obs), 
	                        [data.M[i, j] for (i, j) in data.ids_obs])

	# keep track of performance metric.
	perf_metrics = zeros(nfolds, ngrid)
	fill!(perf_metrics, NaN)

	fig_losses = Figure()
	axs = [Axis(fig_losses[i, j]) for i = 1:nfolds, j = 1:ngrid]
	# loop thru folds
	for (i_fold, (ids_cv_train, ids_cv_test)) in enumerate(kf_split)
		# get the list of matrix entries, vector of tuples
		cv_train_entries = data.ids_obs[ids_cv_train]
		cv_test_entries  = data.ids_obs[ids_cv_test]
		
		###
		# create copy of data
		# introduce additional missing values, where the cv-test data are.
		cv_data = deepcopy(data)
		for (i, j) in cv_test_entries
			cv_data.M[i, j] = missing
			cv_data.M[j, i] = missing
		end
		cv_data.ids_missing = vcat(data.ids_missing, cv_test_entries)
		cv_data.ids_obs = cv_train_entries
		check_consistency(cv_data)
		
		# loop thru hyperparams
		for (n, hps) in enumerate(cv_hyperparams)
			# train model
			cv_model, cv_losses = construct_train_model(hps, data, nb_epochs)
			#_viz_loss!(axs[i_fold, n], cv_losses)
			# evaluate on cv-test data
			perf_metric = compute_perf_metrics(cv_model, cv_data, cv_test_entries)
			perf_metrics[i_fold, n] = perf_metric[type_of_perf_metric]
			#perf_metric[i_fold, n] = accuracy_score(ms_true, ms_pred)
		end		
	end

	opt_hyperparams = cv_hyperparams[argmax(mean(eachrow(perf_metrics))[:])]

	return perf_metrics, opt_hyperparams, fig_losses
end

# ╔═╡ 53987486-ed22-424e-b744-6033ac4be4c6
cv_hyperparams = [(k=rand([2, 3]), γ=rand(Uniform(0, 0.1)), λ=rand(), σ=rand())
				   for _ = 1:ngrid]

# ╔═╡ 434193a3-3b06-494d-b7db-09f152ad437f
perf_metric, opt_hyperparams, fig_losses = do_hyperparam_optimization(data, cv_hyperparams)

# ╔═╡ f1da061e-a8f2-4074-92cf-6183e50e10ba
md"## viz results

train deployment model, `opt_model`, on all obs data using optimal hyperparams.
analyze it more closely. 

1. plot latent vectors (or PCA of them if in higher dims); color by miscibility, different markers for different classes.
2. confusion matrix on test data
3. plot loss as function of iters

and repeat with and without graph regularization.
"

# ╔═╡ 2ccfb38a-4205-429d-9edc-76a20082d735
opt_model, opt_losses = construct_train_model(opt_hyperparams, data, nb_epochs)

# ╔═╡ 1d7cc2ed-f383-4a0b-85aa-849f3f95f983
class_to_marker = Dict("Polymer"    => :circle, 
					   "Protein"    => :rect, 
					   "Surfactant" => :diamond, 
					   "Salt"       => :hexagon)

# ╔═╡ a120c609-f5cb-4012-9f7c-e3702269b541
class_to_color = Dict(zip(["Polymer", "Protein", "Surfactant", "Salt"], ColorSchemes.Accent_4))

# ╔═╡ d39f41a6-e48e-40b7-8929-f62d8ce22a2f
function viz_latent_space(model::Model)
	do_pca = size(model.C)[1] != 2
	
	if do_pca
		# standardize for PCA
		C̃ = deepcopy(opt_model.C')
		for c = 1:size(C̃)[2]
			C̃[:, c] = (C̃[:, c] .- mean(C̃[:, c])) / std(C̃[:, c])
		end
		# the projection into 2D via PCA
		pca = PCA(n_components=2)
		Ĉ = pca.fit_transform(C̃)'
	else
		Ĉ = model.C
	end
	
	fig = Figure()
	ax  = Axis(fig[1, 1], 
		xlabel=do_pca ? "PC1" : "latent dimension 1", 
		ylabel=do_pca ? "PC2" : "latent dimension 2",
		aspect=DataAspect()
	)
	
	vlines!(0.0, color="gray")
	hlines!(0.0, color="gray")
	for class in ["Polymer", "Protein", "Surfactant", "Salt"]
		ids = compound_classes .== class
		scatter!(Ĉ[1, ids], Ĉ[2, ids], 
			strokewidth=1, strokecolor="black", marker=class_to_marker[class],
			color=class_to_color[class], label=class, aspect=DataAspect()
		)
	end
	axislegend()
	return fig
end

# ╔═╡ b93b5880-3a89-4028-a2dc-b93d5c6b138d
viz_latent_space(opt_model)

# ╔═╡ 80ff40cf-1a99-4035-bc97-5e2f4a85c6b0
# on tests data
function compute_cm(model::Model, data::MiscibilityData)
	m = [M_complete[i, j]            for (i, j) in data.ids_missing]
	m̂ = [pred_mᵢⱼ(model, i, j) > 0.5 for (i, j) in data.ids_missing]

	cm = confusion_matrix(m, m̂)
end

# ╔═╡ 976e29ae-393a-48df-9dc2-e393af5dcc7a
function viz_confusion(cm::Matrix)
	cm_to_plot = reverse(cm, dims=1)'

	fig = Figure()
	ax  = Axis(fig[1, 1], 
		xlabel="prediction", ylabel="truth",
		xticks=(1:2, ["immiscible", "miscible"]),
		yticks=(1:2, reverse(["immiscible", "miscible"]))
	)
	# TODO this is messed up.
	hm = heatmap!(cm_to_plot, 
		colormap=ColorSchemes.algae, 
		colorrange=(0, maximum(cm))
	)
	for i = 1:2
        for j = 1:2
            text!("$(Int(cm_to_plot[i, j]))",
                  position=(i, j), align=(:center, :center), color="white", 
				  fontsize=50
			)
        end
    end
    Colorbar(fig[1, 2], hm, label="# pairs")
	return fig
end

# ╔═╡ 4edc1e98-24b6-4cc0-823e-db7137e9ac9a
cory = 5.0

# ╔═╡ a6e3639c-3974-4556-8f45-7fc69a1cd426
function eraser(x)
	cory = 6.0
end

# ╔═╡ 0fb063c6-ff71-4b7b-9ae0-e4022443f046


# ╔═╡ 35666424-d5dc-4a2f-a650-3241a0952c07
viz_confusion(cm) # TODO compare to random guessing

# ╔═╡ 0fd7342c-a459-40dc-9fd8-c8b1d8103c61
md"plot distribution of predictions on test data"

# ╔═╡ 380153ec-1faa-4359-b558-28119edef2ad


# ╔═╡ 258a1618-997e-4ba6-a315-4a679a3ba22d
function viz_preds_on_test(model)
	m = [M_complete[i, j]            for (i, j) in data.ids_missing]
	m̂ = [pred_mᵢⱼ(model, i, j) for (i, j) in data.ids_missing]
	
	fig = Figure()
	ax  = Axis(fig[1, 1], xlabel="ĉ", ylabel="# test data")
	for l in [0, 1]
		hist!(m̂[m .== l], label="class $l", color=(label_to_color[l], 0.5))
	end
	axislegend(position=:lt)
	ylims!(0, nothing)
	return fig
end

# ╔═╡ d508d3f9-ac1a-463b-9930-aa403df25558
viz_preds_on_test(model)

# ╔═╡ 8d664f82-2436-4f6a-9769-f468a1828d6d
M_complete

# ╔═╡ 69ae7968-e6a7-46e7-8506-f0ea2b30e047
data.ids_missing

# ╔═╡ c39c6e79-b5ec-4bf3-8c82-3ba5d043dc51
md"## big function"

# ╔═╡ 8c7201bb-58d3-43f0-a5eb-cce612f1f1d7
function run_experiment(θ::Float64; 
                        nfolds::Int=3, 
                        ngrid::Int=10, 
					    nb_epochs::Int=1000,
                        graph_regularization::Bool=true,
                        class_kernel::Bool=true,
					    seed::Int=97330
)
	# introduce missing values into the matrix.
	data = sim_data_collection(θ, M_complete, seed=seed)

	# create list of hyperparams to explore.
	Random.seed!(seed)
	cv_hyperparams = [
		(k = rand([2, 3, 4, 5]), 
		 γ = graph_regularization ? rand(Uniform(0, 0.1)) : 0.0, 
		 λ = rand(), 
		 σ = class_kernel ? nothing : rand(Uniform(1e-1, 2))
		)
		    for _ = 1:ngrid]

	# conduct hyper-param optimization viz K-folds cross validation on training data
	perf_metric, opt_hyperparams, fig_losses = do_hyperparam_optimization(
		          data, cv_hyperparams, nfolds=nfolds)

	# train deployment model on all training data with opt hyper-params
	opt_model, opt_losses = construct_train_model(opt_hyperparams, data, nb_epochs)

	# compute performance metrics on train and test
	perf_metrics = (test =compute_perf_metrics(opt_model, data, data.ids_missing),
		            train=compute_perf_metrics(opt_model, data, data.ids_obs))

	# confusion matrix on test
	cm = compute_cm(opt_model, data)

	return (data=data, opt_model=opt_model, perf_metrics=perf_metrics, 
		    opt_hyperparams=opt_hyperparams, fig_losses=fig_losses, cm=cm)
end

# ╔═╡ 0b47f501-95fe-48ca-9627-8db4c764e2e3
results_class_kernel = run_experiment(0.5, class_kernel=true)

# ╔═╡ dd97d9ee-c755-4bdc-86fd-58cebbb638f1
viz_confusion(results_class_kernel.cm)

# ╔═╡ bb2745c2-ff9d-48aa-b1a9-baf455dabbbd
results_features = run_experiment(0.5, class_kernel=false)

# ╔═╡ e339ad29-f74f-4c3c-be0c-e28c67d920a3
md"show confusion matrix and similarity matrix and clustering for a _typical_ example (median F1 score)"

# ╔═╡ b836b38e-b90e-4e7b-99fd-4c8eedf88676
md"bar plot of recall, acc, pre, rec over 10 runs"

# ╔═╡ f58e52b9-5504-4b6e-95ea-2d4019ccb5a5


# ╔═╡ 528cc8ca-63b0-47fc-ac29-92f25e04224a
# ╠═╡ disabled = true
#=╠═╡
data2, opt_model2, perf_metric2, opt_hyperparams2, fig_losses2 = run_experiment(.5, class_kernel = false)
  ╠═╡ =#

# ╔═╡ 347f6868-6db2-4ba7-b672-b996e71b4b5b
#=╠═╡
perf_metric2[argmax(mean(eachcol(perf_metric2))[:])]
  ╠═╡ =#

# ╔═╡ f526f24c-0e69-4cfe-9cc7-cc9ec60bc9f1
#=╠═╡
cm2 = compute_cm(opt_model2, data2)
  ╠═╡ =#

# ╔═╡ 0da70997-abca-4d2a-96d7-9800e75cffbf
#=╠═╡
compute_metric(opt_model2, data2)
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
KernelFunctions = "ec8451be-7e33-11e9-00cf-bbf324bd1392"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
MLJBase = "a7f614a8-145f-11e9-1d2a-a57a1082229d"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
ScikitLearn = "3646fa90-6ef7-5e7e-9f22-8aca16db6324"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
AlgebraOfGraphics = "~0.6.13"
CSV = "~0.10.8"
CairoMakie = "~0.10.0"
ColorSchemes = "~3.20.0"
DataFrames = "~1.4.4"
Distributions = "~0.25.79"
KernelFunctions = "~0.10.51"
LogExpFunctions = "~0.3.19"
MLJBase = "~0.21.3"
PlutoUI = "~0.7.49"
ScikitLearn = "~0.6.5"
StatsBase = "~0.33.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "f04ffbe296aafabf7e2bb56bdfc5a5d72b2e1daa"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "52b3b436f8f73133d7bc3a6c71ee7ed6ab2ab754"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.3"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.AlgebraOfGraphics]]
deps = ["Colors", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "KernelDensity", "Loess", "Makie", "PlotUtils", "PooledArrays", "RelocatableFolders", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "3faa74d8af5a9bdcf79236c2984db71ba9e81e74"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.6.13"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "8c73e96bd6817c2597cfd5615b91fca5deccf1af"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.8"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA", "SnoopPrecompile"]
git-tree-sha1 = "a1889ac0cfd046d62404ac3e0a1cb718575ee017"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.10.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "5084cc1a28976dd1642c9f337b28a3cb03e0f7d2"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.7"

[[deps.CategoricalDistributions]]
deps = ["CategoricalArrays", "Distributions", "Missings", "OrderedCollections", "Random", "ScientificTypes", "UnicodePlots"]
git-tree-sha1 = "23fe4c6668776fedfd3747c545cd0d1a5190eb15"
uuid = "af321ab8-2d2e-40a6-b165-3d674595d28e"
version = "0.1.9"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e08915633fcb3ea83bf9d6126292e5bc5c739922"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.13.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d4f69885afa5e6149d0cab3818491565cf41446d"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "e82c3c97b5b4ec111f3c1b55228cebc7510525a2"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.25"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "a7756d098cbabec6b3ac44f369f74915e8cfd70a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.79"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "9a0472ec2f5409db243160a8b030f94c380167a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "38a92e40157100e796690421e34a11c107205c86"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "993c2b4a9a54496b6d8e265db1244db418f37e01"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.1"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "884477b9886a52a84378275737e2823a5c98e349"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.8.1"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "6872f5ec8fd1a38880f027a26739d42dcda6691f"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.2"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "fb28b5dc239d0174d7297310ef7b84a11804dfab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "fe9aea4ed3ec6afdfbeb5a4f39a2208909b162a6"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.5"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "678d136003ed5bceaab05cf64519e3f956ffa4ba"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "342f789fd041a55166764c351da1710db97ce0e0"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.6"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "0cf92ec945125946352f3d46c96976ab972bde6f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.3.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "82aec7a3dd64f4d9584659dc0b62ef7db2ef3e19"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.2.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "9816b296736292a80b9a3200eb7fbb57aaa3917a"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.5"

[[deps.KernelFunctions]]
deps = ["ChainRulesCore", "Compat", "CompositionsBase", "Distances", "FillArrays", "Functors", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Random", "Requires", "SpecialFunctions", "Statistics", "StatsBase", "TensorCore", "Test", "ZygoteRules"]
git-tree-sha1 = "0e9ed147ebbd8c5dc4a05ed164068beeca77dffb"
uuid = "ec8451be-7e33-11e9-00cf-bbf324bd1392"
version = "0.10.51"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "46efcea75c890e5d820e670516dc156689851722"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.5.4"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LossFunctions]]
deps = ["InteractiveUtils", "Markdown", "RecipesBase"]
git-tree-sha1 = "53cd63a12f06a43eef6f4aafb910ac755c122be7"
uuid = "30fc2ffe-d236-52d8-8643-a9d8f7c094a7"
version = "0.8.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MLJBase]]
deps = ["CategoricalArrays", "CategoricalDistributions", "ComputationalResources", "Dates", "DelimitedFiles", "Distributed", "Distributions", "InteractiveUtils", "InvertedIndices", "LinearAlgebra", "LossFunctions", "MLJModelInterface", "Missings", "OrderedCollections", "Parameters", "PrettyTables", "ProgressMeter", "Random", "ScientificTypes", "Serialization", "StatisticalTraits", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "645ad8980fbd61321dc16dc072f97099d9cf60c9"
uuid = "a7f614a8-145f-11e9-1d2a-a57a1082229d"
version = "0.21.3"

[[deps.MLJModelInterface]]
deps = ["Random", "ScientificTypesBase", "StatisticalTraits"]
git-tree-sha1 = "c8b7e632d6754a5e36c0d94a4b466a5ba3a30128"
uuid = "e80e1ace-859a-464e-9ed9-23947d8ae3ea"
version = "1.8.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "MiniQhull", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Serialization", "Setfield", "Showoff", "SignedDistanceFields", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "7154536d78dcde1c4321b50e0e8dda90995f1f6f"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.0"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "5357b0696f7c245941389995e193c127190d45f8"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.0"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.MarchingCubes]]
deps = ["SnoopPrecompile", "StaticArrays"]
git-tree-sha1 = "ffc66942498a5f0d02b9e7b1b1af0f5873142cdc"
uuid = "299715c1-40a9-479a-aaf9-4a633d36f717"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test", "UnicodeFun"]
git-tree-sha1 = "f04120d9adf4f49be242db0b905bea0be32198d1"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.4"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.MiniQhull]]
deps = ["QhullMiniWrapper_jll"]
git-tree-sha1 = "9dc837d180ee49eeb7c8b77bb1c860452634b0d1"
uuid = "978d7f02-9e05-4691-894f-ae31a51d76ca"
version = "0.4.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "5ae7ca23e13855b3aba94550f26146c01d259267"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "f71d8950b724e9ff6110fc948dff5a329f901d64"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "f809158b27eba0c18c269cf2a2be6ed751d3e81d"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.17"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "1155f6f937fa2b94104162f01fa400e192e4272f"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.4.2"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "84a314e3926ba9ec66ac097e3635e270986b0f10"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.9+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "b64719e8b4504983c7fca6cc9db3ebc8acc2a4d6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "5b7690dd212e026bbab1860016a6601cb077ab66"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "96f6db03ab535bdb901300f88335257b0018689d"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "53b8b07b721b77144a0fbbbc2675222ebf40a02d"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.94.1"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QhullMiniWrapper_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Qhull_jll"]
git-tree-sha1 = "607cf73c03f8a9f83b36db0b86a3a9c14179621f"
uuid = "460c41e3-6112-5d7f-b78c-b6823adb3f2d"
version = "1.0.0+1"

[[deps.Qhull_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "238dd7e2cc577281976b9681702174850f8d4cbc"
uuid = "784f63db-0788-585a-bace-daefebcd302b"
version = "8.0.1001+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "18c35ed630d7229c5584b945641a73ca83fb5213"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
git-tree-sha1 = "bc12e315740f3a36a6db85fa2c0212a848bd239e"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.2"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.ScientificTypes]]
deps = ["CategoricalArrays", "ColorTypes", "Dates", "Distributions", "PrettyTables", "Reexport", "ScientificTypesBase", "StatisticalTraits", "Tables"]
git-tree-sha1 = "75ccd10ca65b939dab03b812994e571bf1e3e1da"
uuid = "321657f4-b219-11e9-178b-2701a2544e81"
version = "3.0.2"

[[deps.ScientificTypesBase]]
git-tree-sha1 = "a8e18eb383b5ecf1b5e6fc237eb39255044fd92b"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "3.0.0"

[[deps.ScikitLearn]]
deps = ["Compat", "Conda", "DataFrames", "Distributed", "IterTools", "LinearAlgebra", "MacroTools", "Parameters", "Printf", "PyCall", "Random", "ScikitLearnBase", "SparseArrays", "StatsBase", "VersionParsing"]
git-tree-sha1 = "de6a32950c170e5fd5a2d8bcba0fb97a1028ab06"
uuid = "3646fa90-6ef7-5e7e-9f22-8aca16db6324"
version = "0.6.5"

[[deps.ScikitLearnBase]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "7877e55c1523a4b336b433da39c8e8c08d2f221f"
uuid = "6e75b9c4-186b-50bd-896f-2d2496a4843e"
version = "0.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "efd23b378ea5f2db53a55ae53d3133de4e080aa9"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.16"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "ffc098086f35909741f71ce21d03dadf0d2bfa76"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.11"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "30b9236691858e13f167ce829490a68e1a597782"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "3.2.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "a5e15f27abd2692ccb61a99e0854dfb7d48017db"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.33"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "b03a3b745aa49b566f128977a7dd1be8711c5e71"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.14"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "f8cd5b95aae14d3d88da725414bdde342457366f"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.2"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "e4bdc63f5c6d62e80eb1c0043fcc0360d5950ff7"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.10"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.UnicodePlots]]
deps = ["ColorSchemes", "ColorTypes", "Contour", "Crayons", "Dates", "FileIO", "FreeType", "LinearAlgebra", "MarchingCubes", "NaNMath", "Printf", "Requires", "SnoopPrecompile", "SparseArrays", "StaticArrays", "StatsBase", "Unitful"]
git-tree-sha1 = "e20b01d50cd162593cfd9691628c830769f68987"
uuid = "b8865327-cd53-5732-bb35-84acbb429228"
version = "3.3.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d670a70dd3cdbe1c1186f2f17c9a68a7ec24838c"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.12.2"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╠═c3f6fd50-1697-11ed-1af7-c717d141dfac
# ╠═57504276-6927-496e-9441-79ed5fbdd78e
# ╠═2cbce39d-e9c3-474a-9a0f-ceb6442c27de
# ╠═7b533d68-4bb8-4b1d-a48f-0adffc99ae97
# ╠═b683bdc4-9956-4176-9c03-bebf15d0243e
# ╠═57a6f607-cd48-4e57-a452-06e0fcdd435c
# ╠═9cf628c4-3194-423a-90fc-94bb46c25450
# ╠═5df72253-8573-46e9-a21b-923e0ce3f360
# ╠═cabd0924-cae4-4992-bf49-92353aeea388
# ╠═adb50317-7042-4894-8f2a-e9f6505adcfe
# ╟─f924129b-4e85-45b4-95ce-6846d69d16f3
# ╠═b9e723cf-3af9-40b7-a9de-61893c633d53
# ╠═e559a0c8-7e1f-4fd6-a029-35b6354cd609
# ╠═61377868-d14c-4927-a6e1-6d887903fbf8
# ╠═7b9c60c4-b7a7-499b-8888-508d3610a966
# ╟─35c79c1f-f086-419f-a2dc-2c1710e3ba59
# ╠═dc065d67-229c-4f3f-8912-45a2deb3bf5e
# ╠═22531e96-f0f7-45c0-a327-53094a83ddcf
# ╠═01890165-4cb1-4917-959d-9a0c28e60073
# ╠═9b0388ac-a640-4ffa-8531-b1ae7b68393f
# ╠═8deb8f6f-2be0-4c48-bf2b-4bf2e391c5c1
# ╠═acf912f1-51cf-4943-8964-c51cc3bf3427
# ╠═3a6f1f39-b598-4ce2-8858-b17551133bfe
# ╟─3153dd82-8be6-4683-b93d-9629c3d3312f
# ╠═d6d83535-fa87-4de9-8ddc-0920a6bad075
# ╠═e2e9ae6a-f156-4a91-a63c-6ab7365fe19b
# ╠═e6608357-deeb-4bc4-b5b6-5b5b50f07a3d
# ╠═4aea5499-5cdc-41e4-ad84-ae0eaf6c69e5
# ╠═32837121-174a-498a-9586-a5e4ac13398c
# ╠═1be7e102-d7e9-4c08-989b-211ae302754e
# ╟─ee675570-94ef-4061-9438-7e60a0e5dae1
# ╠═d1062c18-909e-49f7-b7b8-82a566316b64
# ╠═6bc8a921-999d-4be4-a35c-aa3a0383c53e
# ╠═885c20b3-6376-4ce3-991a-8e1db6e88771
# ╠═35e3ff6c-7e08-4956-b6d2-e9c6b8b076c6
# ╠═4c217dfa-2f9c-4ac9-8dfa-3f6bf270139f
# ╟─54b331da-d6c7-4b5d-b71a-bc7d33c3c71a
# ╠═c106b08f-2025-4531-90e2-ed1d4d8e8b36
# ╟─eaa4560e-cb0a-47bc-87ba-9ff4e3faa2c2
# ╠═fabe9de6-f0e3-49fd-84b5-fe973e39f5a0
# ╠═af50a2c4-d8aa-49bd-b596-39964d2d8c7e
# ╠═500a3554-0bbf-42d2-ae4b-0b6c59dbeea6
# ╠═4758b6e9-3c01-42b4-85bf-30ac2298e511
# ╠═81f11413-bc17-4d38-880b-1afbca4ee01d
# ╠═c44a5395-8143-4e33-9eff-dd7f01a1217b
# ╟─a87cc2ef-b91f-49a2-9fbe-cef50192976d
# ╠═8b9046d8-8a0b-4f4a-b461-8c20349ecd30
# ╠═52deabd5-bc18-4659-a7e5-9f42b1e5f591
# ╠═8db228a7-bf2b-4f1e-ab7e-596b2957498f
# ╠═660f9ba7-871a-43df-b69d-d792a8e75456
# ╠═ba7a3de9-b37f-4c67-869b-fceb7915ffab
# ╠═1a178491-815e-49e3-b423-832a24abdba9
# ╠═3f1084c4-28f8-4f0d-98b8-2a0426c89e18
# ╠═0ede5258-7b6a-4975-8d6b-65639bb7ac74
# ╟─c939cb70-858d-4b57-977a-693097c78499
# ╠═6ceaf279-5e77-4585-812b-506c152fd6ab
# ╠═92cd854b-fb84-4981-ba8e-b24ecd1509c7
# ╠═53349731-ff98-4ed7-9877-fe8cf389e6e3
# ╠═5efedaad-892f-4005-8ea2-27ce5b3fd7d8
# ╠═9eef7d3f-a121-41ff-828b-042910b56789
# ╠═53987486-ed22-424e-b744-6033ac4be4c6
# ╠═434193a3-3b06-494d-b7db-09f152ad437f
# ╟─f1da061e-a8f2-4074-92cf-6183e50e10ba
# ╠═2ccfb38a-4205-429d-9edc-76a20082d735
# ╠═1d7cc2ed-f383-4a0b-85aa-849f3f95f983
# ╠═a120c609-f5cb-4012-9f7c-e3702269b541
# ╠═d39f41a6-e48e-40b7-8929-f62d8ce22a2f
# ╠═b93b5880-3a89-4028-a2dc-b93d5c6b138d
# ╠═80ff40cf-1a99-4035-bc97-5e2f4a85c6b0
# ╠═976e29ae-393a-48df-9dc2-e393af5dcc7a
# ╠═4edc1e98-24b6-4cc0-823e-db7137e9ac9a
# ╠═a6e3639c-3974-4556-8f45-7fc69a1cd426
# ╠═0fb063c6-ff71-4b7b-9ae0-e4022443f046
# ╠═35666424-d5dc-4a2f-a650-3241a0952c07
# ╟─0fd7342c-a459-40dc-9fd8-c8b1d8103c61
# ╠═380153ec-1faa-4359-b558-28119edef2ad
# ╠═258a1618-997e-4ba6-a315-4a679a3ba22d
# ╠═d508d3f9-ac1a-463b-9930-aa403df25558
# ╠═8d664f82-2436-4f6a-9769-f468a1828d6d
# ╠═69ae7968-e6a7-46e7-8506-f0ea2b30e047
# ╟─c39c6e79-b5ec-4bf3-8c82-3ba5d043dc51
# ╠═8c7201bb-58d3-43f0-a5eb-cce612f1f1d7
# ╠═0b47f501-95fe-48ca-9627-8db4c764e2e3
# ╠═dd97d9ee-c755-4bdc-86fd-58cebbb638f1
# ╠═bb2745c2-ff9d-48aa-b1a9-baf455dabbbd
# ╟─e339ad29-f74f-4c3c-be0c-e28c67d920a3
# ╟─b836b38e-b90e-4e7b-99fd-4c8eedf88676
# ╠═f58e52b9-5504-4b6e-95ea-2d4019ccb5a5
# ╠═528cc8ca-63b0-47fc-ac29-92f25e04224a
# ╠═347f6868-6db2-4ba7-b672-b996e71b4b5b
# ╠═f526f24c-0e69-4cfe-9cc7-cc9ec60bc9f1
# ╠═0da70997-abca-4d2a-96d7-9800e75cffbf
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
