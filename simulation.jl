### A Pluto.jl notebook ###
# v0.19.20

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
	
	using Revise, PlutoUI, Distributions, DataFrames, ProgressMeter, CairoMakie
	using MiscibilityMF
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

# ╔═╡ 59986921-6e8c-4d38-ae3c-7e7ba1bc246d
md"# read in raw data"

# ╔═╡ 11ed0c10-24d7-4c4e-a9f3-4625d96ee7ff
raw_data = retreive_raw_data()

# ╔═╡ 2090c5f8-cd9c-4bd5-abd5-b5724420c940
raw_data.X

# ╔═╡ 417910f2-15e6-4db5-b7f6-b14275ef43c7
unique(raw_data.classes)

# ╔═╡ 6a683231-7767-4383-a4f0-531b0724513f
raw_data.classes

# ╔═╡ c4a26802-2a30-4b8e-b057-2fce989f0aa0
sum(raw_data.M_complete .== 0) / 2 # number immiscible

# ╔═╡ e6a9d962-21c3-4823-85f4-1dc11a19a251
sum(raw_data.M_complete .== 1) / 2 - raw_data.n_compounds # number miscible

# ╔═╡ f1c5c51c-1d5a-496f-b000-0dc0c8f3383d
viz_miscibility_matrix(raw_data.M_complete, raw_data)

# ╔═╡ cdf7421a-aed6-47fb-aac7-74136355a0d3
md"# introduce missing values"

# ╔═╡ 6a1a696c-88e5-46b3-abcc-376ec8099d90
θ = 0.4 # fraction missing

# ╔═╡ 3c44a682-8161-4b03-aaf0-4d9b813c99cb
data = sim_data_collection(θ, raw_data, weigh_classes=false)

# ╔═╡ e292d20e-9031-4276-87fb-b59685db2f72
length(data.ids_obs)

# ╔═╡ 6991d2b0-73df-4b09-aadf-a5877fa5aa5a
length(data.ids_missing)

# ╔═╡ 8d1b193b-a66c-4f37-80ae-7083f93ceb78
viz_miscibility_matrix(data.M, raw_data)

# ╔═╡ d2913b8d-ccef-4790-aa69-56106767f592
md"# dev model

to handle imbalanced, classes: do not weigh classes, rather adjust threshold for balanced accuracy. 

## hyper-param search
"

# ╔═╡ 5b7bcc88-3048-4701-9349-6de5db12be92
nb_epochs = 250

# ╔═╡ 024ce284-90f3-43e6-b701-1f13d209462f
hyperparams_cv = [(
	k=rand([2, 3]), 
	γ=rand(Uniform(0, 0.1)),
	λ=rand(),
	σ=nothing, 
	use_features=false)
			   for _ = 1:5]

# ╔═╡ 4c53ea02-bd91-44a7-9a32-d4759021b7f8
perf_metrics, opt_hyperparams, fig_losses = do_hyperparam_optimization(data, hyperparams_cv, raw_data, nb_epochs=nb_epochs)

# ╔═╡ 09165856-9117-47e4-8560-fc3f457ad6df
fig_losses

# ╔═╡ 925791d7-3dee-4d6b-9baa-9ee85afb487c
md"## train model with opt hyper-params

TODO: study sensitivity w.r.t. loss"

# ╔═╡ a3457343-dc4d-4046-83d3-b7bdc20c427c
model, losses = construct_train_model(opt_hyperparams, data, raw_data, nb_epochs, record_loss=true, α=0.01 / 1.5)

# ╔═╡ 2533b992-e881-417e-ba9b-7c3555751340
losses[end]

# ╔═╡ a58e7958-458f-4900-8de3-f4eeae945710
viz_loss(losses) # TODO ask Joshua

# ╔═╡ 33e459c3-0d6a-4999-816f-54069bbad86f
begin
	set_opt_cutoff!(model, raw_data, data.ids_obs)
	# model.cutoff = 0.5
end

# ╔═╡ 143582f4-83dc-4f38-befb-eb0109c37b7f
viz_latent_space(model, raw_data)

# ╔═╡ 1069ec41-4733-4111-becd-043a104d1c35
md"

class 1 = miscible = positive (arbitrarily). AND the majority class.

recall = P(pred + | +)

precision = P(+ | pred +)

does not even consider the negative, P(test - | -) = 1 - P(test + | -)
"

# ╔═╡ cf67313d-8b14-40e1-abfd-ab559450e098
perf = compute_perf_metrics(model, raw_data, data.ids_missing)

# ╔═╡ 29400905-c956-4272-b721-9896a41ffce1
perf.cm

# ╔═╡ 77b9fb17-63e6-4c2c-b0ee-e5919358b4dc
viz_confusion(perf.cm, save_fig=true)

# ╔═╡ 1ce32f38-22e1-43a1-8aa8-49141573208c
md"## baseline model"

# ╔═╡ ea729838-bdd8-4e17-823d-ac027de562c8
baseline = test_perf_baseline_model(data, raw_data, set_opt_cutoff=true)

# ╔═╡ f78eccb9-a59f-4d05-a477-bd5ca7ae2e24
viz_confusion(Float64.(baseline.cm) / 2)

# ╔═╡ 2a85c371-731f-4110-b50a-3d196184f8bb
md"# multiple runs and sparsities"

# ╔═╡ 05bd9cc6-ce9e-4dec-8c82-d7c62d4d0b6f
@bind do_multiple_runs CheckBox(default=false)

# ╔═╡ e8221855-745c-4eb1-8320-d785b89c284f
θs = [0.2, 0.5, 0.8]

# ╔═╡ ea48a8dd-d504-4025-bab4-b2f57e1fd256
if do_multiple_runs
	θ_to_perf = Dict()
	@showprogress for θ in θs
		θ_to_perf[θ] = run_experiments(θ, raw_data, nruns=10, nb_hyperparams=15)
	end
end

# ╔═╡ 5906e98f-2d7d-4416-abf0-7d64e927bb40
if do_multiple_runs
	for θ in θs
		# performance
		p = Gadfly.plot(θ_to_perf[θ], x=:metric, y=:score, color=:model,
			Gadfly.Geom.boxplot,
			Gadfly.Guide.title("θ = $θ"),
			Gadfly.Coord.cartesian(ymin=0.4, ymax=1.0)
		)
		p |> Gadfly.PDF("results_$θ.pdf")
	end
end

# ╔═╡ d9ae5b52-22c1-4559-9084-ec38bf3fb4a7
function viz_hyperparam(hp::Symbol, perf_data::DataFrame)# performance
	return Gadfly.plot(filter(row -> row["metric"] == "f1", perf_data), y=hp, x=:model, color=:model,
		Gadfly.Geom.beeswarm
		# Gadfly.Coord.cartesian(ymin=0.4, ymax=1.0)
	)
end

# ╔═╡ 3c8a8663-d588-44b8-a244-e84e7ef964dc
viz_hyperparam(:λ, θ_to_perf[0.2])

# ╔═╡ c09d65a7-3457-4a5a-8521-d01513797dd2
viz_hyperparam(:γ, θ_to_perf[0.2])

# ╔═╡ f37c6347-17fe-4166-bcad-3f5951c70dcb
viz_hyperparam(:k, θ_to_perf[0.2])

# ╔═╡ Cell order:
# ╠═4305ab70-e080-11ed-1f7c-1b8fb559b6c3
# ╠═006c7e92-e8a5-471e-92be-4d27134295f4
# ╠═818f126f-61b2-4a8f-8b47-7b7dae87ba47
# ╠═78ea2361-f2b7-4e0b-ad73-f4ddffecdff7
# ╟─59986921-6e8c-4d38-ae3c-7e7ba1bc246d
# ╠═11ed0c10-24d7-4c4e-a9f3-4625d96ee7ff
# ╠═2090c5f8-cd9c-4bd5-abd5-b5724420c940
# ╠═417910f2-15e6-4db5-b7f6-b14275ef43c7
# ╠═6a683231-7767-4383-a4f0-531b0724513f
# ╠═c4a26802-2a30-4b8e-b057-2fce989f0aa0
# ╠═e6a9d962-21c3-4823-85f4-1dc11a19a251
# ╠═f1c5c51c-1d5a-496f-b000-0dc0c8f3383d
# ╟─cdf7421a-aed6-47fb-aac7-74136355a0d3
# ╠═6a1a696c-88e5-46b3-abcc-376ec8099d90
# ╠═3c44a682-8161-4b03-aaf0-4d9b813c99cb
# ╠═e292d20e-9031-4276-87fb-b59685db2f72
# ╠═6991d2b0-73df-4b09-aadf-a5877fa5aa5a
# ╠═8d1b193b-a66c-4f37-80ae-7083f93ceb78
# ╟─d2913b8d-ccef-4790-aa69-56106767f592
# ╠═5b7bcc88-3048-4701-9349-6de5db12be92
# ╠═024ce284-90f3-43e6-b701-1f13d209462f
# ╠═4c53ea02-bd91-44a7-9a32-d4759021b7f8
# ╠═09165856-9117-47e4-8560-fc3f457ad6df
# ╟─925791d7-3dee-4d6b-9baa-9ee85afb487c
# ╠═a3457343-dc4d-4046-83d3-b7bdc20c427c
# ╠═2533b992-e881-417e-ba9b-7c3555751340
# ╠═a58e7958-458f-4900-8de3-f4eeae945710
# ╠═33e459c3-0d6a-4999-816f-54069bbad86f
# ╠═143582f4-83dc-4f38-befb-eb0109c37b7f
# ╟─1069ec41-4733-4111-becd-043a104d1c35
# ╠═cf67313d-8b14-40e1-abfd-ab559450e098
# ╠═29400905-c956-4272-b721-9896a41ffce1
# ╠═77b9fb17-63e6-4c2c-b0ee-e5919358b4dc
# ╟─1ce32f38-22e1-43a1-8aa8-49141573208c
# ╠═ea729838-bdd8-4e17-823d-ac027de562c8
# ╠═f78eccb9-a59f-4d05-a477-bd5ca7ae2e24
# ╟─2a85c371-731f-4110-b50a-3d196184f8bb
# ╠═05bd9cc6-ce9e-4dec-8c82-d7c62d4d0b6f
# ╠═e8221855-745c-4eb1-8320-d785b89c284f
# ╠═ea48a8dd-d504-4025-bab4-b2f57e1fd256
# ╠═5906e98f-2d7d-4416-abf0-7d64e927bb40
# ╠═d9ae5b52-22c1-4559-9084-ec38bf3fb4a7
# ╠═3c8a8663-d588-44b8-a244-e84e7ef964dc
# ╠═c09d65a7-3457-4a5a-8521-d01513797dd2
# ╠═f37c6347-17fe-4166-bcad-3f5951c70dcb
