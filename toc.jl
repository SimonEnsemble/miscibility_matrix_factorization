### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 3501318c-0035-11ee-1fd1-ffd79845d836
begin
	import Pkg; Pkg.activate()
	push!(LOAD_PATH, "src/")
	
	using CairoMakie, ColorSchemes, Random, LinearAlgebra, Printf, Graphs, GraphMakie, NetworkLayout
end

# ╔═╡ a6f849fe-3a37-48ce-bb1c-6d902ea49747
begin
	import AlgebraOfGraphics
	AlgebraOfGraphics.set_aog_theme!(fonts=[AlgebraOfGraphics.firasans("Light"), AlgebraOfGraphics.firasans("Light")])
	the_resolution = (500, 380)
	update_theme!(
		fontsize=25, 
		linewidth=4,
		markersize=14,
		titlefont=AlgebraOfGraphics.firasans("Light"),
		resolution=the_resolution
	)
end

# ╔═╡ f2abfddc-2a83-4c12-b071-84aefff6d665
begin
	n = 10
	M   = ones(n, n)
	M_c = similar(M)
	for i = 1:n
		for j = i+1:n
			M[i, j] = M[j,  i] = M_c[i, j] = M_c[j, i] = rand([0, 1])
			if rand() < 0.4
				M[i, j] = M[j,  i] = NaN
			end
		end
	end
	M
end

# ╔═╡ e4828911-6c24-4dd8-92b7-a4a7ec65b487
M

# ╔═╡ 5259f765-625f-43bd-95cf-aded752e5558
M_c

# ╔═╡ 7827f454-70fc-447c-98da-5331d3ee350d
classes = vcat(["polymer" for _ = 1:3],
				["surfactant" for _ = 4:8],
				["salt" for i = 9:10])

# ╔═╡ 6b4cac48-bcb5-4535-a6a1-6fafd4bc6984
class_to_color = Dict(zip(unique(classes), ColorSchemes.Paired_5[2:4]))

# ╔═╡ df62dc83-d6b8-4e62-819b-5c25d2354333
function draw_matrix(M)
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

	for i = 1:n
		for j = 1:n
			if isnan(reverse(M', dims=2)[i, j])
				text!(ax, i, j, text="?", align=(:center, :center), font=AlgebraOfGraphics.firasans("Light"))
			end
		end
	end

	@assert length(classes) == n
	c0 = 0.5
	if any(isnan.(M))
		for c in unique(classes) # loop thru classes
			
			l = sum(classes .== c) # number of instances of this class
	
			# draw brackets on top
			bracket!(t_ax, c0 + 0.5, 0, c0 + l - 0.5, 0, orientation=:up, 
				fontsize=17,#fontsize=big_fontsize/1.8,
					 font=AlgebraOfGraphics.firasans("Light"), text=lowercase(c), color=class_to_color[c])
	
			# draw brackets on right
			bracket!(r_ax, 0, n + 1 - (c0 + l) + 0.5, 0, n + 1 - c0 - 0.5,
				fontsize=17,
					 orientation=:down,# fontsize=big_fontsize/1.8,
					 font=AlgebraOfGraphics.firasans("Light"), text=lowercase(c), color=class_to_color[c])
	
			c0 += l
		end

		linkxaxes!(ax, t_ax)
		linkyaxes!(ax, r_ax)
		ylims!(t_ax, 0, 2) # to see text
		xlims!(r_ax, 0, 2) # to see text
	end
	colsize!(fig.layout, 1, Aspect(1, 1.0))

	
	## legend
	legend_patches = [
		PolyElement(color=miscibility_colormap[1], strokecolor="gray", polystrokewidth=1),
		PolyElement(color=miscibility_colormap[end], strokecolor="gray", polystrokewidth=1)
	]
	legend_labels = ["immiscible", "miscible"]
	legend_patchsize = (20, 20)
	# if any(ismissing.(M))
		# push!(legend_patches, PolyElement(color="white", strokecolor="gray", polystrokewidth=1))
		# push!(legend_labels, "missing")
		# legend_patchsize = (14, 14, 14)
	# end
	resize_to_layout!(fig)
#     # this messes up the brackets
	Legend(fig[2, 1], legend_patches, legend_labels, patchsize=legend_patchsize, orientation=:horizontal)
	if any(isnan.(M))
		notify(t_ax.finallimits)
		notify(r_ax.finallimits)
	end
	if any(isnan.(M))
		save("toc_M.pdf", fig)
	else
		save("toc_M_complete.pdf", fig)
	end
	fig
end

# ╔═╡ bd41e64a-6805-4611-bf84-ce1b33747244
draw_matrix(M)

# ╔═╡ 4248935c-9f4c-4923-bf14-6135b60c78af
draw_matrix(M_c)

# ╔═╡ 33cd9c06-ae76-49c5-9942-23d4a9cfec42
begin
	graph = SimpleGraph(n)
	for i = 1:n
		for j = 1:n
			if i == j
				continue
			end
			if classes[i] == classes[j]
				add_edge!(graph, i, j)
			end
		end
	end

	f, ax, p = graphplot(graph, node_size=30,
		node_color=[class_to_color[c] for c in classes], layout=Spring(C=10)
		)
	hidedecorations!(ax); hidespines!(ax)
	save("toc_G.pdf", f)
	f
end

# ╔═╡ 53db750c-52c3-48af-acc8-b24babe63964
classes

# ╔═╡ Cell order:
# ╠═3501318c-0035-11ee-1fd1-ffd79845d836
# ╠═a6f849fe-3a37-48ce-bb1c-6d902ea49747
# ╠═f2abfddc-2a83-4c12-b071-84aefff6d665
# ╠═e4828911-6c24-4dd8-92b7-a4a7ec65b487
# ╠═5259f765-625f-43bd-95cf-aded752e5558
# ╠═7827f454-70fc-447c-98da-5331d3ee350d
# ╠═6b4cac48-bcb5-4535-a6a1-6fafd4bc6984
# ╠═df62dc83-d6b8-4e62-819b-5c25d2354333
# ╠═bd41e64a-6805-4611-bf84-ce1b33747244
# ╠═4248935c-9f4c-4923-bf14-6135b60c78af
# ╠═33cd9c06-ae76-49c5-9942-23d4a9cfec42
# ╠═53db750c-52c3-48af-acc8-b24babe63964
