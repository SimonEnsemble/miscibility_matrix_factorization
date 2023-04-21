"""
build compound miscibility matrix.
* `0` ⟹ immisible
* `1` ⟹  miscible.

[link to original paper](https://pubs.acs.org/doi/10.1021/acsami.0c21036).

`pairs.csv` is \"Supplemental Table 4 006: Sample Code (XLSX)\"

`compounds.csv` is from Supp Table 1 link (but improperly labeled as Supp Table 2...)
"""
function build_miscibility_matrix()
    # miscibility data
    df = CSV.read(joinpath("data", "pairs.csv"), DataFrame, header=0)

    M_complete = Int.(Array(df)) # miscibility matrix
    n_compounds = size(M_complete)[1]
    
    # checks
    @assert M_complete' == M_complete   # symmetric
	@assert all(diag(M_complete) .== 1) # diagonal all ones
	# miscible entries
	@assert M_complete[2, 1] == M_complete[4, 1] == M_complete[13, 3] == 1
	# immiscible entries
	@assert M_complete[4, 3] == M_complete[13, 4] == M_complete[12, 3] == M_complete[15, 1] == 0
	@assert all(diag(M_complete) .== 1) # all compounds miscible with themselves.
    
    return M_complete, n_compounds
end

"""
retreive compound names, classes, and feature matrix.
"""
function build_compound_info()
    # compound info
    compounds = CSV.read("compounds.csv", DataFrame)
    compound_names = String.(compounds[:, "NAME"])
    compound_classes = String.(compounds[:, "CLASS"])

    # compound features
    feature_names = ["monomer_mw", "XlogP3", "h-bond_donors",
					 "h-bond_acceptors", "complexity", "concentration",
					 "polymer_mw"]
	feature_matrix = Matrix(compounds[:, feature_names])
	# normalize
	for j = 1:size(feature_matrix)[2]
		μ = mean(feature_matrix[:, j])
		σ = std(feature_matrix[:, j])
		feature_matrix[:, j] .= (feature_matrix[:, j] .- μ) / σ
	end

    return compound_names, compound_classes, feature_names, collect(feature_matrix')
end

struct RawData
    M_complete::Matrix{Int64}
    n_compounds::Int
    compounds::Vector{String}
    classes::Vector{String}
    features::Vector{String}
    X::Matrix{Float64}
end

function retreive_raw_data()
    M_complete, n_compounds = build_miscibility_matrix()
    compounds, classes, features, X = build_compound_info()

    # sort acc to class
    ids = sortperm(classes)
    return RawData(M_complete[ids, ids], n_compounds, compounds[ids], classes[ids], features, X[:, ids])
end
