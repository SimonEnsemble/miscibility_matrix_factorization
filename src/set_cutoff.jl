function set_opt_cutoff!(model::MFModel, raw_data::RawData, ids::Vector{Tuple{Int, Int}}; type_of_perf_metric::Symbol=:balanced_accuracy)
	cutoffs = range(0.1, 0.9, length=20)
	scores = zeros(20)
	for (n, cutoff) in enumerate(cutoffs) 
		model.cutoff = cutoff
        scores[n] = compute_perf_metrics(model, raw_data, ids)[type_of_perf_metric]
	end

	model.cutoff = cutoffs[argmax(scores)]
end
