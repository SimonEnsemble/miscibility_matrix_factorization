function class_similarity_matrix(raw_data::RawData)
    n_compounds = raw_data.n_compounds
    K = zeros(n_compounds, n_compounds)
    for i = 1:n_compounds
        for j = 1:n_compounds
            if i == j
                K[i, j] = NaN # for safety
            else
                K[i, j] = raw_data.classes[i] == raw_data.classes[j]
            end
        end
    end
    return K
end

function kernel_matrix(raw_data::RawData, σ::Float64)
    kernel = with_lengthscale(SqExponentialKernel(), σ)
    n_compounds = raw_data.n_compounds
    @assert size(raw_data.X)[2] == n_compounds

    K = zeros(n_compounds, n_compounds)
    for i = 1:n_compounds
        for j = 1:n_compounds
            if i == j
                K[i, j] = NaN # for safety
            else
                K[i, j] = kernel(raw_data.X[:, i], raw_data.X[:, j])
            end
        end
    end
    return K
end
