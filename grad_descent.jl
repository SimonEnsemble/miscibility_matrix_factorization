### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ c3f6fd50-1697-11ed-1af7-c717d141dfac
using CSV, DataFrames,Statistics, LogExpFunctions, PyPlot,LinearAlgebra, Random , StatsBase

# ╔═╡ 43e0bde9-7910-46c7-a4d1-e3ec84eeb26e
using PyCall

# ╔═╡ f924129b-4e85-45b4-95ce-6846d69d16f3
md"## read in miscibility data"

# ╔═╡ b9e723cf-3af9-40b7-a9de-61893c633d53
begin
	df = CSV.read("pairs.csv", DataFrame, header=0)
	M_complete = Array(df) # miscibility matrix
end

# ╔═╡ dc065d67-229c-4f3f-8912-45a2deb3bf5e
@assert M_complete' == M_complete # assert symmetric

# ╔═╡ a5e6d4ac-d1c4-44d5-9ea8-a599893376de
@assert all(diag(M_complete) .== 1) # assert diagonal all ones

# ╔═╡ 964f0898-3bb8-4cc4-9fb7-8e073e09c7a5


# ╔═╡ b656b1bf-ab58-4831-929d-95961a7b1505
n_compounds = size(M_complete)[1] #number of compounds

# ╔═╡ 22531e96-f0f7-45c0-a327-53094a83ddcf
begin
	compounds = CSV.read("compounds.csv", DataFrame)
	# add ATPS (Aqueous Two-Phase System) forming fractions data
	compounds[!,:ATPS_frac] = 1 .- sum.(eachrow(df))./n_compounds
	compounds
end

# ╔═╡ 915c174f-dbd1-42bb-92f0-1ccf7d319258
compound_names = String.(CSV.read("compounds.csv", DataFrame)[:, :NAME])

# ╔═╡ d6d83535-fa87-4de9-8ddc-0920a6bad075
compound_classes = String.(CSV.read("compounds.csv", DataFrame)[:, :CLASS])

# ╔═╡ ee675570-94ef-4061-9438-7e60a0e5dae1
md"## introducing missing values"

# ╔═╡ 885c20b3-6376-4ce3-991a-8e1db6e88771
# θ:fraction of missing entries
# important note: we only pay attn to upper triangle of matrix.
function sim_data_collection(θ::Float64, M_complete::Matrix{Int64})
	n_compounds = size(M_complete)[1]

	# make list of entries in upper diagonal of the matrix
	all_entries = [(i, j) for i = 1:n_compounds for j = i+1:n_compounds]
	@assert length(all_entries) == (n_compounds ^ 2 - n_compounds) ÷ 2

	# number of observed entries
	nb_observed = floor(Int, (1-θ) * length(all_entries))

	# sample observed tuples
    ids_obs = StatsBase.sample(all_entries, nb_observed, replace=false)
	# the rest are unobserved
	ids_unobs = [(i, j) for i = 1:n_compounds for j = i+1:n_compounds 
		             if !((i, j) in ids_obs)]
	
	# construct the matrix with missing values
    M = zeros(Union{Int64, Missing}, n_compounds, n_compounds)
	fill!(M, missing) # default missing
	# fill in observed values
	for (i, j) in ids_obs
		M[i, j] = M_complete[i, j]
		M[j, i] = M_complete[i, j]
		
	end
	
    return M, ids_obs, ids_unobs, all_entries
end

# ╔═╡ 35e3ff6c-7e08-4956-b6d2-e9c6b8b076c6
M, ids_obs, ids_unobs, all_entries = sim_data_collection(0.5, M_complete)

# ╔═╡ 0aa5411e-896f-4487-93c3-9711f71bdb3a
ids_obs

# ╔═╡ eaa4560e-cb0a-47bc-87ba-9ff4e3faa2c2
md"## the loss function and its gradient"

# ╔═╡ 5b66b823-d2c7-4e13-80e9-384e19e149e3
# indicator function
function same_class(i::Int, j::Int, compound_classes::Vector{String})
	return compound_classes[i] == compound_classes[j]
end

# ╔═╡ 45d36237-a7d2-4833-add3-bbc95fc427e1
same_class(1, 1, compound_classes)

# ╔═╡ 500a3554-0bbf-42d2-ae4b-0b6c59dbeea6
function loss(C::Matrix{Float64}, 
	          M::Matrix{Union{Missing, Int64}}, 
	          #ids_obs::Vector{Tuple{Int64, Int64}},
	          γ::Float64)
	l = 0.0
	n_compounds = size(M)[1]
	

	# loop over observed data
	# loop over all pairs
	for c = 1:n_compounds
		js_obs = .! ismissing.(M[c, :])
		x = C[c,:]
		for j=1:n_compounds
			y = C[j,:]
			if js_obs[j]
				l +=  -(M[c,j]*log(logistic(sum(x.*y)))+(1-M[c,j])*log(1-logistic(sum(x.*y)))) + γ*same_class(c, j, compound_classes)*sum((x-y).*transpose(x-y))
			end
		end
	end
	
	return l
end

# ╔═╡ 4758b6e9-3c01-42b4-85bf-30ac2298e511
function grad_loss(C::Matrix{Float64}, 
	           c::Int,
	           M::Matrix{Union{Missing, Int64}}, 
	           #ids_obs::Vector{Tuple{Int64, Int64}},
	           γ::Float64)
	g = 0
	n = size(M)[1]
	x = C[c,:]
	js_obs = .! ismissing.(M[c, :])
	
	for j=1:n
		if js_obs[j]
			y = C[j,:]
			g = g .-((M[c,j]*(1-logistic(sum(x.*y)))-(1-M[c,j])*logistic(sum(x.*y)))).*y .+ 2*γ*same_class(c, j, compound_classes).*(x-y)
		end
	end
	
	return g
end

# ╔═╡ f8f73c2f-e332-4b16-85cf-2729780a13ad
# ╠═╡ disabled = true
#=╠═╡
function grad_descent(C0,M,γ)
	C = C0
	n = size(M)[1]
	for t=1:70
		alpha = .01
		for i=1:n
			C[i,:] = C[i,:] .- alpha.*grad_loss(C,c,M,γ) 
		end
	end
	return C
end
  ╠═╡ =#

# ╔═╡ c44a5395-8143-4e33-9eff-dd7f01a1217b
function grad_descent_step(C,M,γ)
	#C = C0
	n = size(M)[1]
	
	alpha = .007
	for c=1:n
		C[c,:] = C[c,:] .- alpha.*grad_loss(C,c,M,γ) 
	end
	
	return C
end

# ╔═╡ 38e663bb-0de3-46c4-9e72-761ea22bc9a6
# sample simulation
begin
	steps = [300 300]
	l1 = zeros(steps[1])
	l2 = zeros(steps[2])
	k = 2
	C1 = rand(Float64,(n_compounds,k))
	C_γ = rand(Float64,(n_compounds,k))
	
	γ = [0.0 .3]
	
	for t=1:steps[1]
		global C1 = grad_descent_step(C1,M,γ[1])
		l1[t] = loss(C1,M,γ[1])
	end

	for t=1:steps[2]
		global C_γ = grad_descent_step(C_γ,M,γ[2])
		l2[t] = loss(C_γ,M,γ[2])
		#l1[t] = dot((C1*transpose(C1)),M)
		#l2[t] = dot((C2*transpose(C2)),M)
	end
end

# ╔═╡ f1da061e-a8f2-4074-92cf-6183e50e10ba
md"## viz results
1. plot latent vectors (or PCA of them if in higher dims); color by miscibility, different markers for different classes.
2. confusion matrix on test data
3. plot loss as function of iters

and repeat with and without graph regularization.
"

# ╔═╡ 3f1084c4-28f8-4f0d-98b8-2a0426c89e18
begin
	figure()
	xlabel("step")
	ylabel("loss")
	scatter(1:steps[1],l1, label="w/o graph regularization")
	scatter(1:steps[2],l2, label="w/ graph regularization")
	legend()
	gcf()
end

# ╔═╡ d39f41a6-e48e-40b7-8929-f62d8ce22a2f
begin
	figure(figsize=(9, 9))
	#fig, (ax1, ax2) = plt.subplots(1, 2)

	markers = zeros(n_compounds)

	markers_dict = Dict("Polymer" => "s", "Protein" => "o", "Surfactant" => "^", "Salt" => "P")
	
	#for (i,class) in enumerate(compounds[:, :CLASS])
		#markers[i] = markers_dict[class]
	#end
	
	scatter(C_γ[:, 1], C_γ[:, 2], edgecolor="k", c=compounds[:, :ATPS_frac], marker=markers, cmap="gist_rainbow")


	xlabel("latent dimension 1")
	ylabel("latent dimension 2")
	
	#for (i,label) in enumerate(compounds[:, :CLASS])
    	#PyPlot.annotate(label, (Y[1,i], Y[2,i]), fontsize=6)
	#end

	axvline(x=0.0, color="lightgray", zorder=0)
	axhline(y=0.0, color="lightgray", zorder=0)
	colorbar(label="ATPS (Aqueous Two-Phase System) forming fraction", fraction=0.046, pad=0.04)
	gca().set_aspect("equal", "box")
	tight_layout()
	#savefig("prop_latent_space.pdf", format="pdf")
	gcf()
	
end

# ╔═╡ 25140d66-18cc-4f0b-a31b-1ab4816ed73e
sk = pyimport("sklearn.metrics")

# ╔═╡ 976e29ae-393a-48df-9dc2-e393af5dcc7a
function cm(C, M_complete, ids_unobs,b,t)
	
	m_true = zeros(length(ids_unobs))
	m_pred = zeros(length(ids_unobs))
	
	for (k,(i,j)) in enumerate(ids_unobs)
		if abs((C*transpose(C))[i,j]-b) > t
			m_pred[k] = 1
		else
			m_pred[k] = 0
		end
		m_true[k] = M_complete[i,j]
	end
		
	cm = sk.confusion_matrix(m_true,m_pred)
	#return sk.ConfusionMatrixDisplay(cm).plot()
	return cm
end

# ╔═╡ 514c2f36-8b34-4cae-b787-ec51d711fa34
md"## confusion matrix: w/ graph regularization"

# ╔═╡ b52d0cb0-6ebb-4f5d-9762-ade7bc305746
begin
	#fig, (ax0, ax1) = subplots(2, 1)
	b = 0.3
	t = 0.5
	cm1=cm(C_γ, M_complete, ids_unobs,b,t)
	fig1=sk.ConfusionMatrixDisplay(cm1).plot()
	gcf()
	
	#fig = cm(C_γ, M_complete, ids_unobs)
	
	
end

# ╔═╡ c7cdb6a2-9544-480c-9d4f-8024e2ea5678
md"## confusion matrix: w/o graph regularization"

# ╔═╡ b20e6084-afa4-4864-b32a-57614bdde09f
begin
	cm2=cm(C1, M_complete, ids_unobs,b ,t)
	fig2=sk.ConfusionMatrixDisplay(cm2).plot()
	gcf()
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CSV = "~0.10.4"
DataFrames = "~0.21.8"
LogExpFunctions = "~0.3.17"
PyCall = "~1.94.1"
PyPlot = "~2.10.0"
StatsBase = "~0.33.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "JSON", "Missings", "Printf", "Statistics", "StructTypes", "Unicode"]
git-tree-sha1 = "2ac27f59196a68070e132b25713f9a5bbc5fa0d2"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.8.3"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "80ca332f6dcb2508adba68f22f551adb2d00a624"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.3"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "9be8be1d8a6f44b96482c8af52238ea7987da3e3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.45.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataFrames]]
deps = ["CategoricalArrays", "Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "Missings", "PooledArrays", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ecd850f3d2b815431104252575e7307256121548"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "0.21.8"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "d19f9edd8c34760dca2de2b503f969d8700ed288"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "361c2b088575b07946508f135ac556751240091c"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.17"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f8c673ccc215eb50fcadb285f522420e29e69e1c"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "0.4.5"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PooledArrays]]
deps = ["DataAPI"]
git-tree-sha1 = "b1333d4eced1826e15adbdf01a4ecaccca9d353c"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "53b8b07b721b77144a0fbbbc2675222ebf40a02d"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.94.1"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "14c1b795b9d764e1784713941e787e1384268103"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.10.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
deps = ["Pkg"]
git-tree-sha1 = "7b1d07f411bc8ddb7977ec7f377b97b158514fe0"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "0.2.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "db8481cf5d6278a121184809e9eb1628943c7704"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.13"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures", "Random", "Test"]
git-tree-sha1 = "03f5898c9959f8115e30bc7226ada7d0df554ddd"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "0.3.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "d24a825a95a6d98c385001212dc9020d609f2d4f"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.8.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═c3f6fd50-1697-11ed-1af7-c717d141dfac
# ╟─f924129b-4e85-45b4-95ce-6846d69d16f3
# ╠═b9e723cf-3af9-40b7-a9de-61893c633d53
# ╠═dc065d67-229c-4f3f-8912-45a2deb3bf5e
# ╠═a5e6d4ac-d1c4-44d5-9ea8-a599893376de
# ╠═964f0898-3bb8-4cc4-9fb7-8e073e09c7a5
# ╠═b656b1bf-ab58-4831-929d-95961a7b1505
# ╠═22531e96-f0f7-45c0-a327-53094a83ddcf
# ╠═915c174f-dbd1-42bb-92f0-1ccf7d319258
# ╠═d6d83535-fa87-4de9-8ddc-0920a6bad075
# ╟─ee675570-94ef-4061-9438-7e60a0e5dae1
# ╠═885c20b3-6376-4ce3-991a-8e1db6e88771
# ╠═35e3ff6c-7e08-4956-b6d2-e9c6b8b076c6
# ╠═0aa5411e-896f-4487-93c3-9711f71bdb3a
# ╟─eaa4560e-cb0a-47bc-87ba-9ff4e3faa2c2
# ╠═5b66b823-d2c7-4e13-80e9-384e19e149e3
# ╠═45d36237-a7d2-4833-add3-bbc95fc427e1
# ╠═500a3554-0bbf-42d2-ae4b-0b6c59dbeea6
# ╠═4758b6e9-3c01-42b4-85bf-30ac2298e511
# ╠═f8f73c2f-e332-4b16-85cf-2729780a13ad
# ╠═c44a5395-8143-4e33-9eff-dd7f01a1217b
# ╠═38e663bb-0de3-46c4-9e72-761ea22bc9a6
# ╟─f1da061e-a8f2-4074-92cf-6183e50e10ba
# ╠═3f1084c4-28f8-4f0d-98b8-2a0426c89e18
# ╠═d39f41a6-e48e-40b7-8929-f62d8ce22a2f
# ╠═43e0bde9-7910-46c7-a4d1-e3ec84eeb26e
# ╠═25140d66-18cc-4f0b-a31b-1ab4816ed73e
# ╠═976e29ae-393a-48df-9dc2-e393af5dcc7a
# ╟─514c2f36-8b34-4cae-b787-ec51d711fa34
# ╠═b52d0cb0-6ebb-4f5d-9762-ade7bc305746
# ╟─c7cdb6a2-9544-480c-9d4f-8024e2ea5678
# ╠═b20e6084-afa4-4864-b32a-57614bdde09f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
