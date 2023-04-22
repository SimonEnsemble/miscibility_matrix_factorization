mutable struct MFModel
	C::Matrix{Float64} # latent vectors in cols
	b::Float64 # bias
	γ::Float64 # regularization param for graph structure
	λ::Float64 # regularization param on latent reps
    use_features::Bool
	σ::Union{Float64, Nothing} # scale param for kernel
	cutoff::Float64 # classification threshold
end

function pred_mᵢⱼ(model::MFModel, i::Int, j::Int)
	cᵢ = model.C[:, i]
	cⱼ = model.C[:, j]
	return logistic(dot(cᵢ, cⱼ) + model.b)
end

function loss(data::MiscibilityData, model::MFModel, K::Matrix{Float64})
	# code below assumes M all 0's or 1's
	@assert all(skipmissing(data.M) .== 1 .|| skipmissing(data.M) .== 0)
	# latent vectors are columns of C
	@assert size(model.C)[2] == size(data.M)[1] == size(data.M)[2]

	# cross-entropy loss expressing mismatch over observations
	l = 0.0
	for i = 1:data.n_compounds
        # regularization of latent vector
		l += model.λ * sum(model.C[:, i] .^ 2) 
		for j = (i+1):data.n_compounds
			###
			#    graph-regularization term
			#      (looks at ALL pairs)
			###
			cᵢ = model.C[:, i]
			cⱼ = model.C[:, j]
			l += model.γ * K[i, j] * sum((cᵢ - cⱼ) .^ 2)

			###
			#    mis-match term
			#      (looks at pair only if not missing.)
			###
			if ! ismissing(data.M[i, j])
				# model prediction
				m̂ᵢⱼ = pred_mᵢⱼ(model, i, j)

				# truth
				mᵢⱼ = data.M[i, j]

				# increment mismatch loss. wt by class.
				wt = data.class_wt[mᵢⱼ]
				l -= (mᵢⱼ * log(m̂ᵢⱼ) + (1 - mᵢⱼ) * log(1 - m̂ᵢⱼ)) * wt
			end
		end
	end
	return l
end

function ∇_cᵢ(data::MiscibilityData, model::MFModel, c::Int, K::Matrix{Float64})
	∇ = zeros(size(model.C)[1])
	
	∇ += model.λ * 2 * model.C[:, c] # regularize latent rep.
	# cross-entropy loss expressing mismatch over observations
	for j = 1:data.n_compounds
		###
		#   contribution to gradient by reg. term 
		###
		if c != j
			∇ += model.γ * K[c, j] * 2 * (model.C[:, c] - model.C[:, j]) # graph
		end
		
		###
		#   contribution to gradient by mis-matches
		###
		if ! ismissing(data.M[c, j])
			# model prediction
			m̂ᵢⱼ = pred_mᵢⱼ(model, c, j)
			# truth
			mᵢⱼ = data.M[c, j]
			# increment loss
			cⱼ = model.C[:, j]
			∇ += cⱼ * (m̂ᵢⱼ - mᵢⱼ) * data.class_wt[mᵢⱼ]
		end
	end
	
	return ∇
end

function ∇_b(data::MiscibilityData, model::MFModel)
	∇ = 0.0

	# cross-entropy loss expressing mismatch over observations
	for i = 1:data.n_compounds
		for j = (i+1):data.n_compounds
			if ! ismissing(data.M[i, j])
				# model prediction
				m̂ᵢⱼ = pred_mᵢⱼ(model, i, j)
				# truth
				mᵢⱼ = data.M[i, j]
				# increment loss
				∇ += (m̂ᵢⱼ - mᵢⱼ) * data.class_wt[mᵢⱼ]
			end
		end
	end
	return ∇
end

# α: learning rate
function grad_descent_epoch!(data::MiscibilityData, model::MFModel, K::Matrix{Float64}; α::Float64=0.01)
	# update compound vectors
	for c = shuffle(1:data.n_compounds)
		model.C[:, c] -= α * ∇_cᵢ(data, model, c, K)
	end
	# update bias
	model.b -= α * ∇_b(data, model)
	return nothing
end

function fraction_miscible(M::Union{Matrix{Int}, Matrix{Union{Missing, Int}}})
	return sum(skipmissing(M) .== 1) / sum(.! ismissing.(M))
end

function construct_train_model(hyperparams::NamedTuple,
	                           data::MiscibilityData,
                               raw_data::RawData,
	                           nb_epochs::Int;
							   α::Float64=0.01, # learning rate
							   cutoff::Float64=0.5,
                               record_loss::Bool=false)
	# initialize model
	f_miscible = fraction_miscible(data.M)
	b_guess = log(f_miscible / (1 - f_miscible))
	C_guess = 0.5 * rand(Uniform(-1, 1), (hyperparams.k, data.n_compounds))
	model = MFModel(C_guess, b_guess, hyperparams.γ, hyperparams.λ, hyperparams.use_features, hyperparams.σ, cutoff)
    
    if hyperparams.use_features
        K = kernel_matrix(raw_data, hyperparams.σ)
    else
        @assert isnothing(hyperparams.σ)
        K = class_similarity_matrix(raw_data)
    end

	# gradient descent epochs. keep track of loss.
	losses = [NaN for _ = 1:nb_epochs]
	for s = 1:nb_epochs
		grad_descent_epoch!(data, model, K, α=α)
        if record_loss
            losses[s] = loss(data, model, K)
        end
	end
	return model, losses
end
