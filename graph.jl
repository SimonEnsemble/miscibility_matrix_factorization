### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ 13a8cc80-e06d-11ed-3734-dde2e620fda7
begin
	import Pkg; Pkg.activate()
	push!(LOAD_PATH, "src/")

	using MiscibilityMF, Graphs, GraphMakie, NetworkLayout, CairoMakie, Random, DataFrames, Statistics
	import Gadfly
end

# ╔═╡ f72f5754-e3bb-4ac9-a03b-dbbf38ec196d
raw_data = retreive_raw_data()

# ╔═╡ f9d05204-d608-4de8-a9f8-d352e18c592d
function construct_graph(raw_data)
	# construct the graph
	graph = SimpleGraph(raw_data.n_compounds)
	
	# loop over pairs of compounds
	for c = 1:raw_data.n_compounds
		for c′ = (c+1):raw_data.n_compounds
			# connect if compounds are the same class
			if raw_data.classes[c] == raw_data.classes[c′]
				add_edge!(graph, c, c′)
			end
		end
	end

	return graph
end

# ╔═╡ 5710f23f-ff98-4ebb-98f4-3c016ff85d3a
graph = construct_graph(raw_data)

# ╔═╡ e927c34e-0dc7-42cf-a97e-d75d775f46a2
begin
	function custom_layout(graph)
		pin = Dict(
			zip(
				[findfirst(raw_data.classes .== c) for c in unique(raw_data.classes)],
				[(5, -5), (-5, 5), (-5, -5), (5, 5)]
			)
		)
		
		positions = spring(graph, C=7.0, seed=5)#, pin=pin)
		
		xs = [pos[1] for pos in positions]
		ys = [pos[2] for pos in positions]
		ids_polymer = raw_data.classes .== "Polymer"
		α = 3.0
		xs[ids_polymer] *= α
		ys[ids_polymer] *= α
		return Point.(zip(xs, ys))
	end
end

# ╔═╡ cb6669f6-72cb-4ab5-bd4c-303449f9b8c9
raw_data.compounds

# ╔═╡ d95c085a-ebf4-4d57-9371-c94219c4c753
begin
	pin = Dict(
		zip(
			[findfirst(raw_data.classes .== c) for c in unique(raw_data.classes)],
			[nothing, (-5, 5), (-5, -5), (5, 5)]
		)
	)
	layout = Spring(;C=10.0, seed=3, pin=pin)

	# layout = SFDP(; iterations=200, C=20.0, K=0.5, seed=4, pin=pin)
	
	node_color = [MiscibilityMF.class_to_color[class] for class in raw_data.classes]
	nlabels, nlabels_offset, nlabels_align = nonoverlapping_nlabels(layout)
	
	f, ax, p = graphplot(graph, 
		node_color=node_color, 
		edge_width=1, 
		edge_color=("black", 0.5),
		layout=layout, 
		nlabels=nlabels, 
		nlabels_align=nlabels_align,
		nlabels_offset=nlabels_offset,
		node_attr=(; markersize=20, strokewidth=1)
	)
	# hidedecorations!(ax); hidespines!(ax)
	
	f
end

# ╔═╡ dcbf750f-aacc-4a10-b096-f154d8e2c55a
nonoverlapping_nlabels(layout, 0.1)

# ╔═╡ 5e40f523-9b42-450a-8b34-102c6427bacf
nlabels

# ╔═╡ bf799b73-b3ee-4871-a409-ff23eeae8f51
rand([3, 4])

# ╔═╡ 6b71b9c0-97d8-4382-acbc-db2c69aa9377
function draw_graph_class(class::String)
	graph_c, nodelabel_c = construct_graph(raw_data, class)
	
	layout = Spring(;C=10.0, seed=7)
	# layout = SFDP(; iterations=200, C=0.6, K=1.0, seed=4)
	f, ax, p = graphplot(graph_c, 
		node_color=MiscibilityMF.class_to_color[class], 
		edge_width=1, 
		edge_color=("black", 0.5),
		layout=layout, 
		nlabels=nodelabel_c,
		node_attr=(; markersize=20, strokewidth=1))
	# hidedecorations!(ax); hidespines!(ax)
	
	f
end

# ╔═╡ 1fd7b8a0-41a9-4174-bc13-f531e379486e
graph

# ╔═╡ aa0e3378-633d-494a-a86a-13933984dde6
function graph_to_dataframe()
	# pin nodes
	pin = Dict(
		zip(
			[findfirst(raw_data.classes .== c) for c in unique(raw_data.classes)],
			[nothing, (-5, 5), (-5, -5), (5, 5)]
		)
	)
	# compute layout
	layout = Spring(;C=10.0, seed=3)#, pin=pin)
	
	df = DataFrame(
		x=[pos[1] for pos in layout(graph)],
		y=[pos[2] for pos in layout(graph)],
		class=lowercase.(raw_data.classes),
		label=raw_data.compounds
	)

	# stretch out polymers
	polymer_center = mean([
		filter(row -> row["class"] == "polymer", df)[:, c] 
		for c in ["x", "y"]]
	)
	shift_center = [-20, 15]
	α = 7.0
	for i = 1:nrow(df)
		if df[i, "class"] == "polymer"
			# sretch
			df[i, "x"] = shift_center[1] + α * (df[i, "x"] - polymer_center[1])
			df[i, "y"] = shift_center[2] + α * (df[i, "y"] - polymer_center[2])
		end
	end
	return df
end

# ╔═╡ 319ca0bb-040d-46a1-8126-6d87c8b78be0
layout(graph)

# ╔═╡ aa4f2212-15fc-4733-b8f5-84a9ea2b7de9
df = graph_to_dataframe()

# ╔═╡ 001f2366-f988-4a0c-a0b4-7692af705c14
begin
	edge_plots = []
	for class in unique(df[:, "class"])
		df_c = filter(row -> row["class"] == class, df)
		n_c = nrow(df_c)
		
		xs =    [df_c[i, "x"] for i = 1:n_c, j = 1:n_c if i < j]
		ys =    [df_c[i, "y"] for i = 1:n_c, j = 1:n_c if i < j]
		xends = [df_c[j, "x"] for i = 1:n_c, j = 1:n_c if i < j]
		yends = [df_c[j, "y"] for i = 1:n_c, j = 1:n_c if i < j]
		
		push!(edge_plots, 
			Gadfly.layer(x=xs, y=ys, xend=xends, yend=yends, Gadfly.Geom.segment, Gadfly.Theme(default_color=colorant"lightgray"))
		)
	end
end

# ╔═╡ a83a153e-a00b-42fc-8afa-092443deb4ba
edge_plots

# ╔═╡ b630bb04-511f-4f80-8268-e71bd1c99f6e
gp = Gadfly.plot(
	Gadfly.layer(df, x="x", y="y", label="label", color="class",
		Gadfly.Geom.point, Gadfly.Geom.label),
	Gadfly.Guide.xlabel(""), Gadfly.Guide.ylabel(""),
	edge_plots...,
	Gadfly.Theme(background_color=colorant"white"),
	Gadfly.Guide.xticks(ticks=nothing), Gadfly.Guide.yticks(ticks=nothing)
)

# ╔═╡ d6fb82d8-8423-4b59-a2f0-36b389f8ea67
Gadfly.draw(Gadfly.SVG("graph.svg", 6Gadfly.inch, 6Gadfly.inch), gp)

# ╔═╡ Cell order:
# ╠═13a8cc80-e06d-11ed-3734-dde2e620fda7
# ╠═f72f5754-e3bb-4ac9-a03b-dbbf38ec196d
# ╠═f9d05204-d608-4de8-a9f8-d352e18c592d
# ╠═5710f23f-ff98-4ebb-98f4-3c016ff85d3a
# ╠═e927c34e-0dc7-42cf-a97e-d75d775f46a2
# ╠═cb6669f6-72cb-4ab5-bd4c-303449f9b8c9
# ╠═dcbf750f-aacc-4a10-b096-f154d8e2c55a
# ╠═d95c085a-ebf4-4d57-9371-c94219c4c753
# ╠═5e40f523-9b42-450a-8b34-102c6427bacf
# ╠═bf799b73-b3ee-4871-a409-ff23eeae8f51
# ╠═6b71b9c0-97d8-4382-acbc-db2c69aa9377
# ╠═1fd7b8a0-41a9-4174-bc13-f531e379486e
# ╠═aa0e3378-633d-494a-a86a-13933984dde6
# ╠═319ca0bb-040d-46a1-8126-6d87c8b78be0
# ╠═aa4f2212-15fc-4733-b8f5-84a9ea2b7de9
# ╠═001f2366-f988-4a0c-a0b4-7692af705c14
# ╠═a83a153e-a00b-42fc-8afa-092443deb4ba
# ╠═b630bb04-511f-4f80-8268-e71bd1c99f6e
# ╠═d6fb82d8-8423-4b59-a2f0-36b389f8ea67
