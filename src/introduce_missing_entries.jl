mutable struct MiscibilityData
	M::Matrix{Union{Int64, Missing}}
	ids_obs::Vector{Tuple{Int, Int}}
	ids_missing::Vector{Tuple{Int, Int}}
    class_wt::Dict{Int, Float64}
	θ::Float64
end

"""
θ:fraction of missing entries
stratified split.
important note: we only pay attn to upper triangle of matrix.
"""
function sim_data_collection(θ::Float64, raw_data::RawData;
                             seed::Union{Nothing, Int}=nothing)
	if ! isnothing(seed)
		Random.seed!(seed)
	end

	@assert (θ < 1.0) && (θ > 0.0)
    n_compounds = raw_data.n_compounds

	# make list of entries in upper diagonal of the matrix. does not include diagonal.
	ids_all_pairs = [(i, j) for i = 1:n_compounds for j = i+1:n_compounds]
	@assert length(ids_all_pairs) == (n_compounds ^ 2 - n_compounds) ÷ 2

	# number of observed entries
	nb_observed = floor(Int, (1-θ) * length(ids_all_pairs))

	# sample observed tuples. stratified split.
	ms_all_pairs = [raw_data.M_complete[i, j] for (i, j) in ids_all_pairs] # for stratifying
    ids_missing, ids_obs = partition(ids_all_pairs, θ, shuffle=true,
		                             stratify=ms_all_pairs)

	# construct the matrix with missing values
    M = zeros(Union{Int64, Missing}, n_compounds, n_compounds)
	fill!(M, missing) # default missing
	# fill in observed values
	for (i, j) in ids_obs
		M[i, j] = raw_data.M_complete[i, j]
		M[j, i] = raw_data.M_complete[i, j]
	end
	# refill diagonal as miscible
	for i = 1:n_compounds
		M[i, i] = 1
	end
	@assert length(ids_obs) * 2 + n_compounds == sum(.! ismissing.(M))

	# to down-weigh the majority class
	fraction_miscible = mean([M[i, j] for (i, j) in ids_obs])
    class_wt = Dict(0 => 1.0, 1 => (1 - fraction_miscible) / fraction_miscible)

    # checks
	@assert all(ismissing.([M[i, j] for (i, j) in ids_missing])) # ids missing actually missing
	@assert all(ismissing.([M[j, i] for (i, j) in ids_missing]))
    
	@assert all(.! ismissing.([M[i, j] for (i, j) in ids_obs])) # ids obs actually obs
	@assert all(.! ismissing.([M[j, i] for (i, j) in ids_obs]))

    @assert all(skipmissing(M') .== skipmissing(M)) # symmetry
	for i = 1:n_compounds, j = 1:n_compounds
		if ! ismissing(M[i, j])
			@assert M[i, j] == raw_data.M_complete[i, j]
		end
	end

    return MiscibilityData(M, ids_obs, ids_missing, class_wt, θ)
end
