### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 3501318c-0035-11ee-1fd1-ffd79845d836
begin
	import Pkg; Pkg.activate()
	push!(LOAD_PATH, "src/")
	
	using CairoMakie, ColorSchemes, Random, LinearAlgebra, Printf, Graphs, GraphMakie
end

# ╔═╡ a6f849fe-3a37-48ce-bb1c-6d902ea49747
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

# ╔═╡ f2abfddc-2a83-4c12-b071-84aefff6d665
begin
	n = 5
	M = ones(n, n)
	for i = 1:n
		for j = i+1:n
			M[i, j] = M[j,  i] = rand([0, 1])
			if rand() < 1/2
				M[i, j] = M[j,  i] = NaN
			end
		end
	end
	M
end

# ╔═╡ df62dc83-d6b8-4e62-819b-5c25d2354333
function draw_matrix()
	miscibility_colormap = reverse(ColorSchemes.:Egypt)
	 
	fig = Figure(resolution=(400, 400))
	ax  = Axis(fig[1, 1],
		xlabel="solution", ylabel="solution",
		xgridvisible=false,
		ygridvisible=false,
		xticks=([], []), 
		yticks=([], []), 
		xticklabelrotation=π/2
	)
	ax.aspect = DataAspect()
	heatmap!(ax, reverse(M', dims=2), colormap=miscibility_colormap)
	for i = 1:n+1
		hlines!(ax, i - 0.5, color="gray", linewidth=1)
		vlines!(ax, i - 0.5, color="gray", linewidth=1)
	end
	
	t_ax = Axis(fig[1, 1, Top()], height=50, xautolimitmargin=(0, 0))
	r_ax = Axis(fig[1, 1, Right()], width=50, yautolimitmargin=(0, 0))
	for a in [t_ax, r_ax]
		hidedecorations!(a)
		hidespines!(a)
	end

	
	classes = ["category X", "category X", "category X", "category Y", "category Y"]
	@assert length(classes) == n
	class_to_color = Dict(zip(unique(classes), ColorSchemes.Egypt[[2, 3]]))
	c0 = 0.5
	for c in unique(classes) # loop thru classes
		
		l = sum(classes .== c) # number of instances of this class

		# draw brackets on top
		bracket!(t_ax, c0 + 0.5, 0, c0 + l - 0.5, 0, orientation=:up, #fontsize=big_fontsize/1.8,
				 font=AlgebraOfGraphics.firasans("Light"), text=lowercase(c), color=class_to_color[c])

		# draw brackets on right
		bracket!(r_ax, 0, n + 1 - (c0 + l) + 0.5, 0, n + 1 - c0 - 0.5,
				 orientation=:down,# fontsize=big_fontsize/1.8,
				 font=AlgebraOfGraphics.firasans("Light"), text=lowercase(c), color=class_to_color[c])

		c0 += l
	end

	linkxaxes!(ax, t_ax)
	linkyaxes!(ax, r_ax)
	ylims!(t_ax, 0, 2) # to see text
	xlims!(r_ax, 0, 2) # to see text
	colsize!(fig.layout, 1, Aspect(1, 1.0))

	
	## legend
	legend_patches = [
		PolyElement(color=miscibility_colormap[1], strokecolor="gray", polystrokewidth=1),
		PolyElement(color=miscibility_colormap[end], strokecolor="gray", polystrokewidth=1)
	]
	legend_labels = ["immiscible", "miscible"]
	legend_patchsize = (35, 35)
	# if any(ismissing.(M))
		push!(legend_patches, PolyElement(color="white", strokecolor="gray", polystrokewidth=1))
		push!(legend_labels, "missing")
		legend_patchsize = (14, 14, 14)
	# end
	resize_to_layout!(fig)
#     # this messes up the brackets
	Legend(fig[0, 1], legend_patches, legend_labels, patchsize=legend_patchsize, orientation=:horizontal)
	# if draw_brackets
		notify(t_ax.finallimits)
		notify(r_ax.finallimits)
	# end
	save("toc_M.pdf", fig)
	fig
end

# ╔═╡ bd41e64a-6805-4611-bf84-ce1b33747244
draw_matrix()

# ╔═╡ 33cd9c06-ae76-49c5-9942-23d4a9cfec42
begin
	graph = SimpleGraph(5)
	add_edge!(graph, 1, 2)
	add_edge!(graph, 3, 4)
	add_edge!(graph, 4, 5)
	add_edge!(graph, 3, 5)

	f, ax, p = graphplot(graph, node_size=30, 
		node_color=vcat([ColorSchemes.Egypt[3] for i = 1:2], [ColorSchemes.Egypt[2] for i = 3:5])
		)
	hidedecorations!(ax); hidespines!(ax)
	save("toc_G.pdf", f)
	f
end

# ╔═╡ Cell order:
# ╠═3501318c-0035-11ee-1fd1-ffd79845d836
# ╠═a6f849fe-3a37-48ce-bb1c-6d902ea49747
# ╠═f2abfddc-2a83-4c12-b071-84aefff6d665
# ╠═df62dc83-d6b8-4e62-819b-5c25d2354333
# ╠═bd41e64a-6805-4611-bf84-ce1b33747244
# ╠═33cd9c06-ae76-49c5-9942-23d4a9cfec42
