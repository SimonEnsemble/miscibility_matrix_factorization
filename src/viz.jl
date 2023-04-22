class_to_color = Dict(zip(["Polymer", "Protein", "Surfactant", "Salt"], ColorSchemes.Accent_4))
class_to_marker = Dict("Polymer"    => :circle,
					   "Protein"    => :rect,
					   "Surfactant" => :diamond,
					   "Salt"       => :hexagon)

function viz_miscibility_matrix(M, raw_data::RawData)
	hm_colormap = reverse(ColorSchemes.turbid)
	
	fig = Figure(resolution=(1200, 1200))
	ax  = Axis(fig[1, 1], 
		xlabel="compound", ylabel="compound", 
		title="miscibility matrix",
		xgridvisible=false, 
		ygridvisible=false, 
		aspect=DataAspect(),
        xticks=(1:raw_data.n_compounds, raw_data.compounds),
        yticks=(1:raw_data.n_compounds, reverse(raw_data.compounds)),
		xticklabelrotation=π/2
	)
	heatmap!(reverse(M', dims=2), colormap=hm_colormap, aspect=DataAspect())
	for i = 1:raw_data.n_compounds
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

_viz_loss!(ax, losses::Vector{Float64}) = lines!(ax, 1:length(losses), losses)

function viz_loss(losses::Vector{Float64})
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="# epochs", ylabel="loss")
	_viz_loss!(ax, losses)
	fig
end

function viz_latent_space(model::MFModel, raw_data::RawData)
	do_pca = size(model.C)[1] != 2

	if do_pca
		# standardize for PCA
		C̃ = deepcopy(model.C')
		for c = 1:size(C̃)[2]
			C̃[:, c] = (C̃[:, c] .- mean(C̃[:, c])) / std(C̃[:, c])
		end
		# the projection into 2D via PCA
		pca = PCA(n_components=2)
		Ĉ = pca.fit_transform(C̃)'
	else
		Ĉ = model.C
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
		ids = raw_data.classes .== class
		scatter!(Ĉ[1, ids], Ĉ[2, ids],
			strokewidth=1, strokecolor="black", marker=class_to_marker[class],
			color=class_to_color[class], label=class, aspect=DataAspect()
		)
	end
	axislegend()
	return fig
end
