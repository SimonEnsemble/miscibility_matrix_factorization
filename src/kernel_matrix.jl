function class_similarity_matrix(raw_data::RawData)
    n_compounds = raw_data.n_compounds
	K = zeros(n_compounds, n_compounds)
	for i = 1:n_compounds
		for j = 1:n_compounds
			if i == j
				K[i, j] = NaN # for safety
			else
				K[i, j] = compound_classes[i] == compound_classes[j]
			end
		end
	end
	return K
end

function kernel_matrix(raw_data::RawData, σ::Float64)
    kernel = with_lengthscale(SqExponentialKernel(), σ)
    n_compounds = raw_data.n_compounds

	K = zeros(n_compounds, n_compounds)
	for i = 1:n_compounds
		for j = 1:n_compounds
			if i == j
				K[i, j] = NaN # for safety
			else
				K[i, j] = kernel(feature_matrix[:, i], feature_matrix[:, j])
			end
		end
	end
	return K
end
