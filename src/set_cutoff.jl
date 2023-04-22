function set_opt_cutoff!(model::MFModel, raw_data::RawData, ids::Vector{Tuple{Int, Int}})
	cutoffs = range(0.1, 0.9, length=20)
	f1_scores = zeros(20)
	for (n, cutoff) in enumerate(cutoffs) 
		model.cutoff = cutoff
		f1_scores[n] = compute_perf_metrics(model, raw_data, ids).f1
	end

	model.cutoff = cutoffs[argmax(f1_scores)]
end
