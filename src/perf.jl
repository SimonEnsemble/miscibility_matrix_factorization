function compute_perf_metrics(model::MFModel, 
	                          data::MiscibilityData, 
	                          ids::Vector{Tuple{Int, Int}})
	m = [M_complete[i, j]            for (i, j) in ids]
	m̂ = [pred_mᵢⱼ(model, i, j) > model.cutoff for (i, j) in ids]

	return (f1=f1_score(m, m̂),
		    accuracy=accuracy_score(m, m̂),
	        precision=precision_score(m, m̂),
		    recall=recall_score(m, m̂)
	)
end
