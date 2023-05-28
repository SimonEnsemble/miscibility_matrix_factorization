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

# ╔═╡ 6ca3a2da-9cdb-453f-bb64-6f0941be1175
graph = construct_graph(raw_data)

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
	
	f, ax, p = graphplot(graph, 
		node_color=node_color, 
		edge_width=1, 
		edge_color=("black", 0.5),
		layout=layout, 
		nlabels=raw_data.compounds, 
		# nlabels_align=nlabels_align,
		# nlabels_offset=nlabels_offset,
		node_attr=(; markersize=20, strokewidth=1)
	)
	hidedecorations!(ax); hidespines!(ax)
	
	f
end

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
# ╠═6ca3a2da-9cdb-453f-bb64-6f0941be1175
# ╠═d95c085a-ebf4-4d57-9371-c94219c4c753
# ╠═aa0e3378-633d-494a-a86a-13933984dde6
# ╠═319ca0bb-040d-46a1-8126-6d87c8b78be0
# ╠═aa4f2212-15fc-4733-b8f5-84a9ea2b7de9
# ╠═001f2366-f988-4a0c-a0b4-7692af705c14
# ╠═b630bb04-511f-4f80-8268-e71bd1c99f6e
# ╠═d6fb82d8-8423-4b59-a2f0-36b389f8ea67
