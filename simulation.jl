### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 4305ab70-e080-11ed-1f7c-1b8fb559b6c3
begin
	import Pkg; Pkg.activate()
	push!(LOAD_PATH, "src/")
	
	using Revise, PlutoUI, Distributions, DataFrames, ProgressMeter, CairoMakie, ColorSchemes, Random, LinearAlgebra, Printf, JLD2
	using MiscibilityMF
	import MLJBase: train_test_pairs, StratifiedCV
end

# ╔═╡ 006c7e92-e8a5-471e-92be-4d27134295f4
begin
	import AlgebraOfGraphics
	AlgebraOfGraphics.set_aog_theme!(fonts=[AlgebraOfGraphics.firasans("Light"), AlgebraOfGraphics.firasans("Light")])
	the_resolution = (500, 380)
	update_theme!(
		fontsize=20, 
		linewidth=4,
		markersize=14,
		titlefont=AlgebraOfGraphics.firasans("Light"),
		resolution=the_resolution
	)
end

# ╔═╡ 818f126f-61b2-4a8f-8b47-7b7dae87ba47
import Gadfly

# ╔═╡ 78ea2361-f2b7-4e0b-ad73-f4ddffecdff7
TableOfContents()

# ╔═╡ f1f1fa73-08b4-4d1d-ad5e-f0903519fb8b
my_seed = 97330 # to make results reproducible

# ╔═╡ 59986921-6e8c-4d38-ae3c-7e7ba1bc246d
md"# read in raw data"

# ╔═╡ 11ed0c10-24d7-4c4e-a9f3-4625d96ee7ff
const raw_data = retreive_raw_data(normalize_features=false)

# ╔═╡ 298bc445-01cf-4c25-8494-5259f4fa6684
raw_data.classes

# ╔═╡ 2090c5f8-cd9c-4bd5-abd5-b5724420c940
raw_data.X

# ╔═╡ 417910f2-15e6-4db5-b7f6-b14275ef43c7
unique(raw_data.classes)

# ╔═╡ 6a683231-7767-4383-a4f0-531b0724513f
raw_data.classes

# ╔═╡ 1ae9210a-102d-4d55-8e32-9ad0950e0e7f
for c in unique(raw_data.classes)
	@printf("# %s = %d\n", c, sum(raw_data.classes .== c))
end

# ╔═╡ c4a26802-2a30-4b8e-b057-2fce989f0aa0
sum(raw_data.M_complete .== 0) / 2 # number immiscible

# ╔═╡ e6a9d962-21c3-4823-85f4-1dc11a19a251
(sum(raw_data.M_complete .== 1) - raw_data.n_compounds) / 2 # number miscible (discount diagonal)

# ╔═╡ f602a244-4781-4155-9f8f-af89783bead3
raw_data.n_compounds * (raw_data.n_compounds - 1) / 2 # pairs

# ╔═╡ f1c5c51c-1d5a-496f-b000-0dc0c8f3383d
viz_miscibility_matrix(raw_data.M_complete, raw_data, draw_brackets=true, savename="miscibility_matrix_complete.pdf")

# ╔═╡ 7dc9e6eb-095d-417d-9cff-534b702159d2
raw_data.classes

# ╔═╡ 13eeb4b6-3017-4e6e-a872-21afc10ff93c
viz_category_miscibility(raw_data)

# ╔═╡ 4c7abede-63b3-4ae9-a14c-f371ee5c460a
md"some manual checks with the paper"

# ╔═╡ 49229699-6b42-4fa7-ad75-6ebdd0852eaf
id_chit = findfirst(raw_data.compounds .== "Chit")

# ╔═╡ 6d3eb9ac-0f2d-4ae2-bc0f-42b5a4beec82
@assert sum(raw_data.M_complete[id_chit, :] .== 0) == 22

# ╔═╡ a5ad1ba4-c3d6-4694-94c7-daf11ee8fd79
@assert raw_data.X[raw_data.features .== "XlogP3", id_chit][1] ≈ -2.7

# ╔═╡ 9562eb0c-5c08-4be3-b578-720c0f74bee1
@assert raw_data.X[raw_data.features .== "is_Polymer", id_chit][1] ≈ 1.0

# ╔═╡ 04fb2587-6566-4221-b403-8b2e6d25811d
@assert raw_data.classes[id_chit] == "Polymer"

# ╔═╡ 1d7a5810-a610-469a-9653-d52d1c697431
id_cas = findfirst(raw_data.compounds .== "Cas")

# ╔═╡ 80984df5-aa7a-4e4d-b0e3-40d5721d400e
@assert sum(raw_data.M_complete[id_cas, :] .== 0) == 19

# ╔═╡ e4d68be8-0adc-4023-bfd2-ff358230d159
@assert raw_data.M_complete[id_chit, id_cas] == 0

# ╔═╡ 6a6b0516-649d-43a1-8635-97c2cc463009
id_chol = findfirst(raw_data.compounds .== "Chol")

# ╔═╡ 53daf2c2-6df5-472d-b0db-84c0958a2046
@assert raw_data.M_complete[id_chit, id_chol] == 0

# ╔═╡ 2675e9aa-a1c6-4982-a83f-a9ff5e3368b5
id_ppg = findfirst(raw_data.compounds .== "PPG")

# ╔═╡ 8874b2df-b17b-43bd-b16a-9dc18e4ad0e1
@assert raw_data.M_complete[id_chit, id_ppg] == 1

# ╔═╡ cdf7421a-aed6-47fb-aac7-74136355a0d3
md"# introduce missing values"

# ╔═╡ 6a1a696c-88e5-46b3-abcc-376ec8099d90
θ = 0.4 # fraction missing

# ╔═╡ 3c44a682-8161-4b03-aaf0-4d9b813c99cb
data = sim_data_collection(θ, raw_data, weigh_classes=false, seed=my_seed)

# ╔═╡ 377d754f-c6cb-48c6-8ce4-1fe3335a23a1
mean([data.M[i, j] for (i, j) in data.ids_obs]) # fraction miscible

# ╔═╡ e292d20e-9031-4276-87fb-b59685db2f72
length(data.ids_obs)

# ╔═╡ 6991d2b0-73df-4b09-aadf-a5877fa5aa5a
length(data.ids_missing)

# ╔═╡ 436164ee-01b8-46ad-9e96-5e72368444d4
the_compound_labels = rich.(raw_data.compounds)

# ╔═╡ 8d1b193b-a66c-4f37-80ae-7083f93ceb78
viz_miscibility_matrix(data.M, raw_data, draw_brackets=true, savename="miscibility_matrix.pdf")

# ╔═╡ 86e8f156-dd43-496c-b040-b4772fe0c536
raw_data.classes

# ╔═╡ d3c890fd-b7b5-4dad-b937-8d98f5278114
viz_miscibility_matrix(data.M, raw_data, draw_brackets=true, savename="miscibility_matrix_simple.pdf", show_solute_labels=false, huge_font=true)

# ╔═╡ d2913b8d-ccef-4790-aa69-56106767f592
md"# dev model

to handle imbalanced, classes: adjust threshold for balanced accuracy. 

set learning rate and number of epochs for gradient descent.
"

# ╔═╡ abbb1485-8652-4cc5-b749-ab93db6b64fc
begin
	α = 0.02
	nb_epochs = 250
end

# ╔═╡ ff07a8bf-4fe4-46ca-a929-d64b557903d6
md"## hyper-param search"

# ╔═╡ 024ce284-90f3-43e6-b701-1f13d209462f
begin
	nb_hyperparams = 25
	
	hyperparams_cv = gen_hyperparams(nb_hyperparams, true, seed=my_seed) # 2nd arg is graph reg
end

# ╔═╡ 4c53ea02-bd91-44a7-9a32-d4759021b7f8
begin
	nfolds = 3
	perf_metrics, opt_hps_id, opt_hyperparams, fig_losses = do_hyperparam_optimization(data, hyperparams_cv, raw_data, 
		nb_epochs=nb_epochs, α=α, record_loss=true, use_adam=true, nfolds=nfolds,
		seed=my_seed
	)
end

# ╔═╡ 65af0efe-5c23-4bc8-898a-ef7a5b43e479
opt_hyperparams

# ╔═╡ 70a722a3-b0b5-4e96-b371-345ddb12665d
viz_hyperparams(hyperparams_cv, opt_hyperparams)

# ╔═╡ b8d7ddb6-aa97-4660-a391-3b7b6ffa59f9
mean(perf_metrics[:, opt_hps_id]) # best cross-validation perf

# ╔═╡ 09165856-9117-47e4-8560-fc3f457ad6df
fig_losses

# ╔═╡ 925791d7-3dee-4d6b-9baa-9ee85afb487c
md"## train model with opt hyper-params"

# ╔═╡ a3457343-dc4d-4046-83d3-b7bdc20c427c
model, losses = construct_train_model(opt_hyperparams, data, raw_data, nb_epochs, record_loss=true, α=α, use_adam=true, seed=my_seed)

# ╔═╡ a58e7958-458f-4900-8de3-f4eeae945710
viz_loss(losses, save_fig=true)

# ╔═╡ 479c0086-0474-448c-a2f9-8ae75aadcc80
losses[end]

# ╔═╡ 33e459c3-0d6a-4999-816f-54069bbad86f
set_opt_cutoff!(model, raw_data, data.ids_obs)

# ╔═╡ 17ffbe7f-b496-470f-b92c-6a78f8c551e5
viz_imputations(model, data, raw_data)

# ╔═╡ 22b47934-4838-43bc-8e8c-b628db62af7c
md"### latent space"

# ╔═╡ 143582f4-83dc-4f38-befb-eb0109c37b7f
viz_latent_space(model, raw_data, save_fig=true)

# ╔═╡ b38dcbae-6895-466d-bc49-7a6865326469
viz_latent_space_3d(model, raw_data)

# ╔═╡ 7a31f94e-d779-4e67-8333-0346c9445ed4
viz_C(model, raw_data, savename="C.pdf")

# ╔═╡ 00f52f67-8afe-4814-92b0-f08f1544bed7
viz_C(model, raw_data, savename="C_simple.pdf", minimal_viz=true, draw_brackets=false)

# ╔═╡ 137a96cd-d5d8-40b7-bbb0-ab3db6521e14
md"viz the imputed matrix"

# ╔═╡ 91229965-9831-4f9f-9052-a2c8d774c58d
begin
	M_predicted = pred_M(model) .> model.cutoff
	viz_miscibility_matrix(M_predicted, raw_data, 
			savename="miscibility_matrix_imputed_simple", show_solute_labels=false)
	
	viz_miscibility_matrix(M_predicted, raw_data, 
			savename="miscibility_matrix_imputed")
end

# ╔═╡ 1069ec41-4733-4111-becd-043a104d1c35
md"

class 1 = miscible = positive (arbitrarily). AND the majority class.

recall = P(pred + | +)

precision = P(+ | pred +)

does not even consider the negative, P(test - | -) = 1 - P(test + | -)
"

# ╔═╡ 6500c702-835b-4fd3-a54a-ad35e11c4111
md"### performance"

# ╔═╡ cf67313d-8b14-40e1-abfd-ab559450e098
perf = compute_perf_metrics(model, raw_data, data.ids_missing)

# ╔═╡ 29400905-c956-4272-b721-9896a41ffce1
perf.cm

# ╔═╡ 77b9fb17-63e6-4c2c-b0ee-e5919358b4dc
viz_confusion(perf.cm, save_fig=true)

# ╔═╡ 1ce32f38-22e1-43a1-8aa8-49141573208c
md"## baseline models

(random forest)
"

# ╔═╡ ea729838-bdd8-4e17-823d-ac027de562c8
baseline = test_perf_baseline_model(data, raw_data, set_opt_cutoff=true)

# ╔═╡ f78eccb9-a59f-4d05-a477-bd5ca7ae2e24
viz_confusion(Float64.(baseline.cm) / 2) # divide by two b.c doubled

# ╔═╡ 1592a5c8-100d-44f5-b48c-2ced1fe601f2
md"(random guessing)"

# ╔═╡ 9e0e2df4-e71b-4fe7-ab63-e7c3402756df
baseline_guessing = test_perf_guessing(data, raw_data)

# ╔═╡ 2b00c2ae-ceba-4fcd-8ea5-9dae2013ee2b
viz_confusion(baseline_guessing.cm)

# ╔═╡ 19d0a18c-ad97-4f76-8cb5-fb3b0cc830d3
begin
	# check confusion matrix
	@assert length(data.ids_missing) == sum(baseline_guessing.cm)
	@assert sum([raw_data.M_complete[i, j] for (i, j) in data.ids_missing] .== 1) == sum(baseline_guessing.cm[2, :]) # sum of second row is truly miscible
	@assert sum([raw_data.M_complete[i, j] for (i, j) in data.ids_missing] .== 0) == sum(baseline_guessing.cm[1, :]) # sum of first row is truly immisible
end

# ╔═╡ d02ee133-9ee6-4ce0-84db-fcdfb78d6830
md"ordinary LMF with same hyper-params but γ=0"

# ╔═╡ 12825032-0994-4d7c-90ef-e473596fb42b
lmf_hyperparams = (;k=opt_hyperparams.k, λ=opt_hyperparams.λ, use_features=false, σ=nothing, γ=0.0)

# ╔═╡ d48a8ef3-4d6e-4378-a206-0901f8cc28e0
lmf_model, _ = construct_train_model(lmf_hyperparams, data, raw_data, nb_epochs, record_loss=true, α=α, use_adam=true, seed=my_seed)

# ╔═╡ 638b1fa7-b664-4173-8a88-271e7e9ba809
set_opt_cutoff!(lmf_model, raw_data, data.ids_obs)

# ╔═╡ 2d2ee780-6d9f-4ca8-84f3-b652349b290c
lmf_perf = compute_perf_metrics(lmf_model, raw_data, data.ids_missing)

# ╔═╡ 7bcf3a25-b298-4638-8c7c-7e39c4e808c9
viz_confusion(lmf_perf.cm)

# ╔═╡ 515d942e-b72f-40f1-a288-4d89aed61d56
viz_latent_space(lmf_model, raw_data, append_filename="lmf", save_fig=true)

# ╔═╡ d935b22f-5121-41c1-8408-c71c0e53f130
viz_latent_space_3d(lmf_model, raw_data, append_filename="lmf")

# ╔═╡ f5d6f5c5-92ae-435e-9d4a-17890d49bc64
md"## loss without Adam, for comparison"

# ╔═╡ d4ce1f41-c8a3-4524-836e-e3099abfae92
_, _losses = construct_train_model(opt_hyperparams, data, raw_data, nb_epochs, record_loss=true, α=α, use_adam=false, seed=my_seed)

# ╔═╡ 621c8c22-e093-4d55-8e80-1d8b81646e69
viz_loss(_losses, save_fig=true, append_filename="without_adam")

# ╔═╡ 73e89a22-fded-47f6-a821-847b1ec13f9c
md"## RF feature importance"

# ╔═╡ da34e60e-b628-4af6-9d1d-1a64a50119d3
viz_rf_feature_importance(raw_data, θ, 25)

# ╔═╡ 4ab546f3-6596-40c3-a4ba-051e97649dd7
md"## timing"

# ╔═╡ 154429da-6316-4683-99ab-ded882f84e06
@time construct_train_model(opt_hyperparams, data, raw_data, nb_epochs, record_loss=true, α=α, use_adam=true)

# ╔═╡ 74bb7d26-993a-4047-9546-9933a64bc1af
@time do_hyperparam_optimization(data, hyperparams_cv, raw_data, 
		nb_epochs=nb_epochs, α=α, record_loss=true, use_adam=true, nfolds=nfolds)

# ╔═╡ 2a85c371-731f-4110-b50a-3d196184f8bb
md"# multiple runs and sparsities"

# ╔═╡ 05bd9cc6-ce9e-4dec-8c82-d7c62d4d0b6f
@bind do_multiple_runs CheckBox(default=false)

# ╔═╡ e8221855-745c-4eb1-8320-d785b89c284f
begin
	θs = [0.2, 0.5, 0.8]
	nruns = 10
end

# ╔═╡ ea48a8dd-d504-4025-bab4-b2f57e1fd256
if do_multiple_runs
	θ_to_perf = Dict()
	for θ in θs
		mf_settings = (; α=α, nb_hyperparams=nb_hyperparams, nb_epochs=nb_epochs, use_adam=true)
		θ_to_perf[θ] = run_experiments(θ, raw_data, nruns, mf_settings)
	end
	jldsave("multiple_runs.jld2"; θ_to_perf)
else
	θ_to_perf = load_object("results_$(nruns)_runs.jld2")
end

# ╔═╡ fe624b98-6283-42d6-86db-ef641edc3a68
md"standard deviation of performance"

# ╔═╡ 29eafcac-3054-44e5-89b0-a93fc981b534
for θ in θs
	σ = std(filter(
		row -> (row["model"] == "GR-LMF") && (row["metric"] == "balanced_accuracy"), 
		θ_to_perf[θ])[:, "score"])
	@show (θ, σ)
end

# ╔═╡ 112e75a8-73e0-429f-9afd-5d4636b0e31c
md"save results to file"

# ╔═╡ 5906e98f-2d7d-4416-abf0-7d64e927bb40
if do_multiple_runs
	for θ in θs
		# performance
		p = Gadfly.plot(θ_to_perf[θ], x=:metric, y=:score, color=:model,
			Gadfly.Geom.boxplot,
			Gadfly.Guide.title("θ = $θ"),
			Gadfly.Coord.cartesian(ymin=0.4, ymax=1.0)
		)
		p |> Gadfly.PDF("results_$(θ)_$(nruns)_runs.pdf")
	end
end

# ╔═╡ cddb1b28-e5b2-436a-b4e1-decb2fa2aae0
function balanced_acc_boxplot(θs, θ_to_perf)
	models = ["GR-LMF", "LMF", "RF", "guess"]
	model_to_dodge = Dict(zip(models, 1:length(models)))
	model_to_color = Dict(zip(models, ColorSchemes.Accent_4))

	# panel for each θ.
	fig = Figure(resolution=(600, 400))
	axs = [Axis(fig[1, i], 
				xlabel=i == 2 ? "model" : "",
				xticklabelrotation=π/2,
				ylabel=i == 1 ? "balanced accuracy" : "",
				xticks=(1:length(models), models),
				title="θ = $(θs[i])"
			)
			for i = 1:length(θs)
	]
	linkyaxes!(axs...)
	ylims!(0.39, 1)
	hideydecorations!(axs[2])
	hideydecorations!(axs[3])
	hidespines!(axs[2], :l)
	hidespines!(axs[3], :l)
	for (i, θ) in enumerate(θs)
		# get data, but only balanced acc
		perf_this_theta = deepcopy(θ_to_perf[θ])
		filter!(row -> row["metric"] == "balanced_accuracy", perf_this_theta)
		for model in models
			d = filter(row -> row["model"] == model, perf_this_theta)
			xs = [model_to_dodge[model] for _ = 1:nrow(d)]
			ys = d[:, "score"]
			boxplot!(axs[i], xs, ys, 
				medianlinewidth=2, mediancolor="black", color=model_to_color[model])#, color = map(d -> dodge_to_color[d], dodge))
		end
	end
	save("balanced_acc.pdf", fig)
	return fig
end

# ╔═╡ 43beda5f-7bcc-4fb4-8b4e-d995bb4c7985
if do_multiple_runs
	balanced_acc_boxplot(θs, θ_to_perf)
end

# ╔═╡ d9ae5b52-22c1-4559-9084-ec38bf3fb4a7
function viz_hyperparam(hp::Symbol, perf_data::DataFrame)# performance
	return Gadfly.plot(
		filter(row -> row["metric"] == "f1", perf_data), 
		y=hp, x=:model, color=:model,
		Gadfly.Geom.beeswarm
		# Gadfly.Coord.cartesian(ymin=0.4, ymax=1.0)
	)
end

# ╔═╡ 3c8a8663-d588-44b8-a244-e84e7ef964dc
if do_multiple_runs
	viz_hyperparam(:λ, θ_to_perf[0.8])
end

# ╔═╡ c09d65a7-3457-4a5a-8521-d01513797dd2
if do_multiple_runs
	viz_hyperparam(:γ, θ_to_perf[0.8])
end

# ╔═╡ f37c6347-17fe-4166-bcad-3f5951c70dcb
if do_multiple_runs
	viz_hyperparam(:k, θ_to_perf[0.8])
end

# ╔═╡ 8f5043d8-0a44-4619-8659-2444acaae5d0
md"# conceptualize our ML procedures
viz cross-validation split
"

# ╔═╡ 317b8e4f-1ad8-47f9-81f2-0271f8e704c7
viz_ml_procedure(data, raw_data, model) # for presentations etc.

# ╔═╡ Cell order:
# ╠═4305ab70-e080-11ed-1f7c-1b8fb559b6c3
# ╠═006c7e92-e8a5-471e-92be-4d27134295f4
# ╠═818f126f-61b2-4a8f-8b47-7b7dae87ba47
# ╠═78ea2361-f2b7-4e0b-ad73-f4ddffecdff7
# ╠═f1f1fa73-08b4-4d1d-ad5e-f0903519fb8b
# ╟─59986921-6e8c-4d38-ae3c-7e7ba1bc246d
# ╠═11ed0c10-24d7-4c4e-a9f3-4625d96ee7ff
# ╠═298bc445-01cf-4c25-8494-5259f4fa6684
# ╠═2090c5f8-cd9c-4bd5-abd5-b5724420c940
# ╠═417910f2-15e6-4db5-b7f6-b14275ef43c7
# ╠═6a683231-7767-4383-a4f0-531b0724513f
# ╠═1ae9210a-102d-4d55-8e32-9ad0950e0e7f
# ╠═c4a26802-2a30-4b8e-b057-2fce989f0aa0
# ╠═e6a9d962-21c3-4823-85f4-1dc11a19a251
# ╠═f602a244-4781-4155-9f8f-af89783bead3
# ╠═f1c5c51c-1d5a-496f-b000-0dc0c8f3383d
# ╠═7dc9e6eb-095d-417d-9cff-534b702159d2
# ╠═13eeb4b6-3017-4e6e-a872-21afc10ff93c
# ╟─4c7abede-63b3-4ae9-a14c-f371ee5c460a
# ╠═49229699-6b42-4fa7-ad75-6ebdd0852eaf
# ╠═6d3eb9ac-0f2d-4ae2-bc0f-42b5a4beec82
# ╠═a5ad1ba4-c3d6-4694-94c7-daf11ee8fd79
# ╠═9562eb0c-5c08-4be3-b578-720c0f74bee1
# ╠═04fb2587-6566-4221-b403-8b2e6d25811d
# ╠═1d7a5810-a610-469a-9653-d52d1c697431
# ╠═80984df5-aa7a-4e4d-b0e3-40d5721d400e
# ╠═e4d68be8-0adc-4023-bfd2-ff358230d159
# ╠═6a6b0516-649d-43a1-8635-97c2cc463009
# ╠═53daf2c2-6df5-472d-b0db-84c0958a2046
# ╠═2675e9aa-a1c6-4982-a83f-a9ff5e3368b5
# ╠═8874b2df-b17b-43bd-b16a-9dc18e4ad0e1
# ╟─cdf7421a-aed6-47fb-aac7-74136355a0d3
# ╠═6a1a696c-88e5-46b3-abcc-376ec8099d90
# ╠═3c44a682-8161-4b03-aaf0-4d9b813c99cb
# ╠═377d754f-c6cb-48c6-8ce4-1fe3335a23a1
# ╠═e292d20e-9031-4276-87fb-b59685db2f72
# ╠═6991d2b0-73df-4b09-aadf-a5877fa5aa5a
# ╠═436164ee-01b8-46ad-9e96-5e72368444d4
# ╠═8d1b193b-a66c-4f37-80ae-7083f93ceb78
# ╠═86e8f156-dd43-496c-b040-b4772fe0c536
# ╠═d3c890fd-b7b5-4dad-b937-8d98f5278114
# ╟─d2913b8d-ccef-4790-aa69-56106767f592
# ╠═abbb1485-8652-4cc5-b749-ab93db6b64fc
# ╟─ff07a8bf-4fe4-46ca-a929-d64b557903d6
# ╠═024ce284-90f3-43e6-b701-1f13d209462f
# ╠═4c53ea02-bd91-44a7-9a32-d4759021b7f8
# ╠═65af0efe-5c23-4bc8-898a-ef7a5b43e479
# ╠═70a722a3-b0b5-4e96-b371-345ddb12665d
# ╠═b8d7ddb6-aa97-4660-a391-3b7b6ffa59f9
# ╠═09165856-9117-47e4-8560-fc3f457ad6df
# ╟─925791d7-3dee-4d6b-9baa-9ee85afb487c
# ╠═a3457343-dc4d-4046-83d3-b7bdc20c427c
# ╠═a58e7958-458f-4900-8de3-f4eeae945710
# ╠═479c0086-0474-448c-a2f9-8ae75aadcc80
# ╠═33e459c3-0d6a-4999-816f-54069bbad86f
# ╠═17ffbe7f-b496-470f-b92c-6a78f8c551e5
# ╟─22b47934-4838-43bc-8e8c-b628db62af7c
# ╠═143582f4-83dc-4f38-befb-eb0109c37b7f
# ╠═b38dcbae-6895-466d-bc49-7a6865326469
# ╠═7a31f94e-d779-4e67-8333-0346c9445ed4
# ╠═00f52f67-8afe-4814-92b0-f08f1544bed7
# ╟─137a96cd-d5d8-40b7-bbb0-ab3db6521e14
# ╠═91229965-9831-4f9f-9052-a2c8d774c58d
# ╟─1069ec41-4733-4111-becd-043a104d1c35
# ╟─6500c702-835b-4fd3-a54a-ad35e11c4111
# ╠═cf67313d-8b14-40e1-abfd-ab559450e098
# ╠═29400905-c956-4272-b721-9896a41ffce1
# ╠═77b9fb17-63e6-4c2c-b0ee-e5919358b4dc
# ╟─1ce32f38-22e1-43a1-8aa8-49141573208c
# ╠═ea729838-bdd8-4e17-823d-ac027de562c8
# ╠═f78eccb9-a59f-4d05-a477-bd5ca7ae2e24
# ╟─1592a5c8-100d-44f5-b48c-2ced1fe601f2
# ╠═9e0e2df4-e71b-4fe7-ab63-e7c3402756df
# ╠═2b00c2ae-ceba-4fcd-8ea5-9dae2013ee2b
# ╠═19d0a18c-ad97-4f76-8cb5-fb3b0cc830d3
# ╟─d02ee133-9ee6-4ce0-84db-fcdfb78d6830
# ╠═12825032-0994-4d7c-90ef-e473596fb42b
# ╠═d48a8ef3-4d6e-4378-a206-0901f8cc28e0
# ╠═638b1fa7-b664-4173-8a88-271e7e9ba809
# ╠═2d2ee780-6d9f-4ca8-84f3-b652349b290c
# ╠═7bcf3a25-b298-4638-8c7c-7e39c4e808c9
# ╠═515d942e-b72f-40f1-a288-4d89aed61d56
# ╠═d935b22f-5121-41c1-8408-c71c0e53f130
# ╟─f5d6f5c5-92ae-435e-9d4a-17890d49bc64
# ╠═d4ce1f41-c8a3-4524-836e-e3099abfae92
# ╠═621c8c22-e093-4d55-8e80-1d8b81646e69
# ╟─73e89a22-fded-47f6-a821-847b1ec13f9c
# ╠═da34e60e-b628-4af6-9d1d-1a64a50119d3
# ╟─4ab546f3-6596-40c3-a4ba-051e97649dd7
# ╠═154429da-6316-4683-99ab-ded882f84e06
# ╠═74bb7d26-993a-4047-9546-9933a64bc1af
# ╟─2a85c371-731f-4110-b50a-3d196184f8bb
# ╠═05bd9cc6-ce9e-4dec-8c82-d7c62d4d0b6f
# ╠═e8221855-745c-4eb1-8320-d785b89c284f
# ╠═ea48a8dd-d504-4025-bab4-b2f57e1fd256
# ╟─fe624b98-6283-42d6-86db-ef641edc3a68
# ╠═29eafcac-3054-44e5-89b0-a93fc981b534
# ╟─112e75a8-73e0-429f-9afd-5d4636b0e31c
# ╠═5906e98f-2d7d-4416-abf0-7d64e927bb40
# ╠═cddb1b28-e5b2-436a-b4e1-decb2fa2aae0
# ╠═43beda5f-7bcc-4fb4-8b4e-d995bb4c7985
# ╠═d9ae5b52-22c1-4559-9084-ec38bf3fb4a7
# ╠═3c8a8663-d588-44b8-a244-e84e7ef964dc
# ╠═c09d65a7-3457-4a5a-8521-d01513797dd2
# ╠═f37c6347-17fe-4166-bcad-3f5951c70dcb
# ╟─8f5043d8-0a44-4619-8659-2444acaae5d0
# ╠═317b8e4f-1ad8-47f9-81f2-0271f8e704c7
