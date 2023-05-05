### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 13a8cc80-e06d-11ed-3734-dde2e620fda7
begin
	import Pkg; Pkg.activate()
	push!(LOAD_PATH, "src/")

	using MiscibilityMF, Graphs, GraphPlot, Compose
end

# ╔═╡ f72f5754-e3bb-4ac9-a03b-dbbf38ec196d
raw_data = retreive_raw_data()

# ╔═╡ 7209567b-c832-460d-ac2f-b241b80774c7
c′=2

# ╔═╡ f9d05204-d608-4de8-a9f8-d352e18c592d
begin
	graph = SimpleGraph(raw_data.n_compounds)
	for c = 1:raw_data.n_compounds
		for c′ = (c+1):raw_data.n_compounds
			if raw_data.classes[c] == raw_data.classes[c′]
				add_edge!(graph, c, c′)
			end
		end
	end
end

# ╔═╡ a82b6265-12e2-4a1a-8e44-de08742dbfcb
begin
	nodelabel = 1:nv(graph)
	nodelabel = raw_data.compounds
	nodefillc = [MiscibilityMF.class_to_color[class] for class in raw_data.classes]
	
	locs_x, locs_y = spring_layout(graph; C=12)

	# omit some labels
	
	nodelabel[1:sum(raw_data.classes .== "Polymer")] .= ""
	p = rand(1:sum(raw_data.classes .== "Polymer"),10)
	nodelabel[p] = raw_data.compounds[p]
	
	#gp = gplot(graph, locs_x, locs_y, nodelabel=nodelabel, nodefillc=nodefillc)
	gp = gplot(graph, nodelabel=nodelabel, nodefillc=nodefillc)
	#draw(PDF("graph.pdf", 16cm, 16cm), gp)
	gp
end

# ╔═╡ Cell order:
# ╠═13a8cc80-e06d-11ed-3734-dde2e620fda7
# ╠═f72f5754-e3bb-4ac9-a03b-dbbf38ec196d
# ╠═7209567b-c832-460d-ac2f-b241b80774c7
# ╠═f9d05204-d608-4de8-a9f8-d352e18c592d
# ╠═a82b6265-12e2-4a1a-8e44-de08742dbfcb
