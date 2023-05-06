class_to_color = Dict(zip(["Polymer", "Protein", "Surfactant", "Salt"], ColorSchemes.seaborn_bright[1:4]))
class_to_marker = Dict("Polymer"    => :circle,
					   "Protein"    => :rect,
					   "Surfactant" => :diamond,
					   "Salt"       => :cross)
miscibility_colormap = reverse(ColorSchemes.turbid)

function viz_miscibility_matrix(M, raw_data::RawData)
    big_fontsize = 50

	fig = Figure(resolution=(1400, 1400))
	ax  = Axis(fig[1, 1], 
		xlabel="compound", ylabel="compound", 
		xgridvisible=false, 
        xlabelsize=big_fontsize,
        ylabelsize=big_fontsize,
		ygridvisible=false, 
		#aspect=DataAspect(),
        xticks=(1:raw_data.n_compounds, raw_data.compounds),
        yticks=(1:raw_data.n_compounds, reverse(raw_data.compounds)),
		xticklabelrotation=π/2
	)
	heatmap!(ax, reverse(M', dims=2), colormap=miscibility_colormap)
	for i = 1:raw_data.n_compounds+1
		hlines!(ax, i - 0.5, color="gray", linewidth=1)
		vlines!(ax, i - 0.5, color="gray", linewidth=1)
	end
	
	# add brackets to indicate the class of the compound
    t_ax = Axis(fig[1, 1, Top()], height=50, xautolimitmargin=(0, 0))
    r_ax = Axis(fig[1, 1, Right()], width=50, yautolimitmargin=(0, 0))
    for a in [t_ax, r_ax]
        hidedecorations!(a)
        hidespines!(a)
    end

	c0 = 0.5
    for c in unique(raw_data.classes) # loop thru classes
		l = sum(raw_data.classes .== c) # number of instances of this class

        # draw brackets on bottom
        bracket!(t_ax, c0 + 0.5, 0, c0 + l - 0.5, 0, orientation=:up, fontsize=big_fontsize/2, font=AlgebraOfGraphics.firasans("Light"), text=lowercase(c), color=class_to_color[c])

        # draw brackets on top
        bracket!(r_ax, 0, c0 + 0.5, 0, c0 + l - 0.5, orientation=:down, fontsize=big_fontsize/2, font=AlgebraOfGraphics.firasans("Light"), text=lowercase(c), color=class_to_color[c])
		
        c0 += l
	end
    linkxaxes!(ax, t_ax)
    linkyaxes!(ax, r_ax)
    ylims!(t_ax, 0, 2) # to see text
    xlims!(r_ax, 0, 2) # to see text
	
	## legend
	legend_patches = [
		PolyElement(color=miscibility_colormap[1], strokecolor="gray"), 
        PolyElement(color=miscibility_colormap[end], strokecolor="gray")
	]
	legend_labels = ["immiscible", "miscible"]
	legend_patchsize = (35, 35)
	if any(ismissing.(M))
		push!(legend_patches, PolyElement(color="white", strokecolor="gray"))
		push!(legend_labels, "missing")
		legend_patchsize = (35, 35, 35)
	end
    # this messes up the brackets
	Legend(fig[0, 1], legend_patches, legend_labels, patchsize=legend_patchsize, orientation=:horizontal, labelsize=big_fontsize)
    save("miscibility_matrix.pdf", fig)
	return fig
end

_viz_loss!(ax, losses::Vector{Float64}) = lines!(ax, 1:length(losses), losses)

function viz_loss(losses::Vector{Float64})
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="# epochs", ylabel="loss")
	_viz_loss!(ax, losses)
    save("loss.pdf", fig)
	fig
end

function viz_latent_space(model::MFModel, raw_data::RawData; incl_legend::Bool=true)
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

	vlines!(0.0, color="lightgray", linewidth=1, linestyle=:dash)
	hlines!(0.0, color="lightgray", linewidth=1, linestyle=:dash)
    sps = []
	classes = ["Polymer", "Protein", "Surfactant", "Salt"]
	for class in classes
		ids = raw_data.classes .== class
		push!(sps,
              scatter!(Ĉ[1, ids], Ĉ[2, ids],
                strokewidth=1, strokecolor="black", marker=class_to_marker[class],
                color=class_to_color[class], label=class, aspect=DataAspect())
             )
	end
    if incl_legend
        Legend(fig[0, 1], sps, classes, orientation=:horizontal)
    end
    save("latent_space.pdf", fig)
	return fig
end

function viz_confusion(cm::Matrix; save_fig::Bool=false)
	cm_to_plot = reverse(cm, dims=1)'

	fig = Figure()
	ax  = Axis(fig[1, 1],
		xlabel="prediction", ylabel="truth",
		xticks=(1:2, ["immiscible", "miscible"]),
		yticks=(1:2, reverse(["immiscible", "miscible"])),
        yticklabelrotation=π/2
	)
	hm = heatmap!(cm_to_plot,
		colormap=ColorSchemes.algae,
		colorrange=(0, maximum(cm))
	)
	for i = 1:2
        for j = 1:2
            text!("$(round(Int, cm_to_plot[i, j]))",
                  position=(i, j), align=(:center, :center), color="white",
				  fontsize=50
			)
        end
    end
    Colorbar(fig[1, 2], hm, label="# pairs")
    if save_fig
        save("confusion_matrix.pdf", fig)
    end
	return fig
end
