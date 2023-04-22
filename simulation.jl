### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ 4305ab70-e080-11ed-1f7c-1b8fb559b6c3
begin
	import Pkg; Pkg.activate()
	push!(LOAD_PATH, "src/")
	
	using Revise, PlutoUI, Distributions
	using MiscibilityMF
end

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

# ╔═╡ c4a26802-2a30-4b8e-b057-2fce989f0aa0
sum(raw_data.M_complete .== 0) / 2 # number immiscible

# ╔═╡ e6a9d962-21c3-4823-85f4-1dc11a19a251
sum(raw_data.M_complete .== 1) / 2 - raw_data.n_compounds # number miscible

# ╔═╡ f1c5c51c-1d5a-496f-b000-0dc0c8f3383d
# viz_miscibility_matrix(raw_data.M_complete, raw_data)

# ╔═╡ cdf7421a-aed6-47fb-aac7-74136355a0d3
md"# introduce missing values"

# ╔═╡ 6a1a696c-88e5-46b3-abcc-376ec8099d90
θ = 0.4 # fraction missing

# ╔═╡ 3c44a682-8161-4b03-aaf0-4d9b813c99cb
data = sim_data_collection(θ, raw_data)

# ╔═╡ e292d20e-9031-4276-87fb-b59685db2f72
length(data.ids_obs)

# ╔═╡ 6991d2b0-73df-4b09-aadf-a5877fa5aa5a
length(data.ids_missing)

# ╔═╡ d2913b8d-ccef-4790-aa69-56106767f592
md"# dev model

## hyper-param search
"

# ╔═╡ 5b7bcc88-3048-4701-9349-6de5db12be92
nb_epochs = 500

# ╔═╡ 024ce284-90f3-43e6-b701-1f13d209462f
hyperparams_cv = [(k=rand([2, 3]), γ=rand(Uniform(0, 0.1)), λ=rand(), σ=nothing, use_features=false)
			   for _ = 1:5]

# ╔═╡ 4c53ea02-bd91-44a7-9a32-d4759021b7f8
perf_metrics, opt_hyperparams, fig_losses = do_hyperparam_optimization(data, hyperparams_cv, raw_data)

# ╔═╡ 09165856-9117-47e4-8560-fc3f457ad6df
fig_losses

# ╔═╡ 925791d7-3dee-4d6b-9baa-9ee85afb487c
md"## train model with opt hyper-params"

# ╔═╡ a3457343-dc4d-4046-83d3-b7bdc20c427c
model, losses = construct_train_model(opt_hyperparams, data, raw_data, nb_epochs)

# ╔═╡ a58e7958-458f-4900-8de3-f4eeae945710
viz_loss(losses)

# ╔═╡ 143582f4-83dc-4f38-befb-eb0109c37b7f
viz_latent_space(model, raw_data)

# ╔═╡ cf67313d-8b14-40e1-abfd-ab559450e098
compute_perf_metrics(model, raw_data, data.ids_missing)

# ╔═╡ 1067cb85-a684-4a11-b4f8-173553e203df
cm = compute_cm(model, raw_data, data.ids_missing)

# ╔═╡ 51fc3d64-452b-41de-b2d2-789015fc4bd3
# TODO set cutoffs

# ╔═╡ 77b9fb17-63e6-4c2c-b0ee-e5919358b4dc
viz_confusion(cm)

# ╔═╡ Cell order:
# ╠═4305ab70-e080-11ed-1f7c-1b8fb559b6c3
# ╠═78ea2361-f2b7-4e0b-ad73-f4ddffecdff7
# ╟─59986921-6e8c-4d38-ae3c-7e7ba1bc246d
# ╠═11ed0c10-24d7-4c4e-a9f3-4625d96ee7ff
# ╠═2090c5f8-cd9c-4bd5-abd5-b5724420c940
# ╠═417910f2-15e6-4db5-b7f6-b14275ef43c7
# ╠═c4a26802-2a30-4b8e-b057-2fce989f0aa0
# ╠═e6a9d962-21c3-4823-85f4-1dc11a19a251
# ╠═f1c5c51c-1d5a-496f-b000-0dc0c8f3383d
# ╟─cdf7421a-aed6-47fb-aac7-74136355a0d3
# ╠═6a1a696c-88e5-46b3-abcc-376ec8099d90
# ╠═3c44a682-8161-4b03-aaf0-4d9b813c99cb
# ╠═e292d20e-9031-4276-87fb-b59685db2f72
# ╠═6991d2b0-73df-4b09-aadf-a5877fa5aa5a
# ╟─d2913b8d-ccef-4790-aa69-56106767f592
# ╠═5b7bcc88-3048-4701-9349-6de5db12be92
# ╠═024ce284-90f3-43e6-b701-1f13d209462f
# ╠═4c53ea02-bd91-44a7-9a32-d4759021b7f8
# ╠═09165856-9117-47e4-8560-fc3f457ad6df
# ╟─925791d7-3dee-4d6b-9baa-9ee85afb487c
# ╠═a3457343-dc4d-4046-83d3-b7bdc20c427c
# ╠═a58e7958-458f-4900-8de3-f4eeae945710
# ╠═143582f4-83dc-4f38-befb-eb0109c37b7f
# ╠═cf67313d-8b14-40e1-abfd-ab559450e098
# ╠═1067cb85-a684-4a11-b4f8-173553e203df
# ╠═51fc3d64-452b-41de-b2d2-789015fc4bd3
# ╠═77b9fb17-63e6-4c2c-b0ee-e5919358b4dc
