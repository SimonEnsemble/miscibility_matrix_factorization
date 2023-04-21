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
