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

function pred_M(model::MFModel)
    M̃ = logistic.(model.C' * model.C .+ model.b)
    # some tests
    @assert M̃[3, 2] ≈ pred_mᵢⱼ(model, 3, 2)
    return M̃
end

function loss(data::MiscibilityData, model::MFModel, K::Matrix{Float64})
    # code below assumes M all 0's or 1's
    @assert all(skipmissing(data.M) .== 1 .|| skipmissing(data.M) .== 0)
    # latent vectors are columns of C
    @assert size(model.C)[2] == size(data.M)[1] == size(data.M)[2]

    l_ce = 0.0 # cross entropy term
    l_nr = 0.0 # norm regularization term
    l_gr = 0.0 # graph-regularization term
    for i = 1:data.n_compounds
        # regularization of latent vector
        l_nr += model.λ * sum(model.C[:, i] .^ 2) 
        for j = (i+1):data.n_compounds
            ###
            #    graph-regularization term
            #      (looks at ALL pairs)
            ###
            cᵢ = model.C[:, i]
            cⱼ = model.C[:, j]
            l_gr += model.γ * K[i, j] * sum((cᵢ - cⱼ) .^ 2)

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
                l_ce -= (mᵢⱼ * log(m̂ᵢⱼ) + (1 - mᵢⱼ) * log(1 - m̂ᵢⱼ)) * wt
            end
        end
    end
    return l_ce + l_nr + l_ce
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

mutable struct AdamOptInfo
    α::Float64 # learning rate
    β::Tuple{Float64, Float64}
    ϵ::Float64
    # moments and biases for compound vectors
    m_C::Matrix{Float64}
    v_C::Matrix{Float64}
    # moments and biases for bias 
    m_b::Float64
    v_b::Float64
end

# see https://arxiv.org/pdf/1412.6980.pdf
function adam_grad_descent_epoch!(data::MiscibilityData, model::MFModel, K::Matrix{Float64}, t::Int, adam_info::AdamOptInfo)
    # update compound vectors
    for c = shuffle(1:data.n_compounds)
        # compute graident
        the_∇_cᵢ = ∇_cᵢ(data, model, c, K)
        # update first moment
        adam_info.m_C[:, c] = adam_info.β[1] * adam_info.m_C[:, c] + (1 - adam_info.β[1]) * the_∇_cᵢ
        # update second moment
        adam_info.v_C[:, c] = adam_info.β[2] * adam_info.v_C[:, c] + (1 - adam_info.β[2]) * the_∇_cᵢ .^ 2
        # bias corrected first moment estimate
        m̂_C = adam_info.m_C[:, c] / (1 - adam_info.β[1] ^ t)
        v̂_C = adam_info.v_C[:, c] / (1 - adam_info.β[2] ^ t)
        # finall, update
        model.C[:, c] -= adam_info.α * m̂_C ./ (sqrt.(v̂_C) .+ adam_info.ϵ)
    end

    # update bias
    the_∇_b = ∇_b(data, model)
    adam_info.m_b = adam_info.β[1] * adam_info.m_b + (1 - adam_info.β[1]) * the_∇_b
    adam_info.v_b = adam_info.β[2] * adam_info.v_b + (1 - adam_info.β[2]) * the_∇_b ^ 2
    m̂_b = adam_info.m_b / (1 - adam_info.β[1] ^ t)
    v̂_b = adam_info.v_b / (1 - adam_info.β[2] ^ t)
    model.b -= adam_info.α * m̂_b / (sqrt(v̂_b) + adam_info.ϵ)
    return nothing
end

function fraction_miscible(M::Union{Matrix{Int}, Matrix{Union{Missing, Int}}})
    return sum(skipmissing(M) .== 1) / sum(.! ismissing.(M))
end

function construct_train_model(hyperparams::NamedTuple,
                               data::MiscibilityData,
                               raw_data::RawData,
                               nb_epochs::Int;
                               α::Float64=0.006, # learning rate
                               cutoff::Float64=0.5,
                               record_loss::Bool=false,
                               use_adam::Bool=false,
                               seed::Union{Nothing, Int}=nothing
                               )
    if ! isnothing(seed)
        Random.seed!(seed)
    end
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

    # scale learning rate with fraction missing entries
    θ̃ = (sum(ismissing.(data.M))) / (prod(size(data.M)) - size(data.M)[1])

    # gradient descent epochs. keep track of loss.
    losses = [NaN for _ = 1:nb_epochs]
    if use_adam
        adam_info = AdamOptInfo(α, (0.9, 0.99), 1e-8, 
                                zeros(size(model.C)), zeros(size(model.C)),
                                0.0, 0.0)
    end
    for t = 1:nb_epochs
        if use_adam
            adam_grad_descent_epoch!(data, model, K, t, adam_info)
        else
            grad_descent_epoch!(data, model, K, α=α * θ̃) # scale learning rate
        end
        if record_loss
            losses[t] = loss(data, model, K)
        end
    end
    return model, losses
end
