function do_hyperparam_optimization(
	data::MiscibilityData,
	hyperparams::Vector{<:NamedTuple},
    raw_data::RawData;
	nfolds::Int=3,
	nb_epochs::Int=250,
	type_of_perf_metric::Symbol=:balanced_accuracy,
    record_loss::Bool=false,
    α::Float64=0.006,
    use_adam::Bool=false
)
	# cross-validation split
	kf_split = train_test_pairs(StratifiedCV(nfolds=nfolds, shuffle=true), 
		                    1:length(data.ids_obs), 
	                        [data.M[i, j] for (i, j) in data.ids_obs])

	# keep track of performance metric.
    perf_metrics = zeros(nfolds, length(hyperparams))
	fill!(perf_metrics, NaN)
	
    fig_losses = Figure(resolution=(3000, 500))
    if record_loss
        axs = [Axis(fig_losses[i, j]) for i = 1:nfolds, j = 1:length(hyperparams)]
    end
	# loop thru folds
	for (i_fold, (ids_cv_train, ids_cv_test)) in enumerate(kf_split)
		# get the list of matrix entries, vector of tuples
		cv_train_entries = data.ids_obs[ids_cv_train]
		cv_test_entries  = data.ids_obs[ids_cv_test]
		
		###
		# create copy of data
		# introduce additional missing values, where the cv-test data are.
		cv_data = deepcopy(data)
		for (i, j) in cv_test_entries
			cv_data.M[i, j] = missing
			cv_data.M[j, i] = missing
		end
		cv_data.ids_missing = vcat(data.ids_missing, cv_test_entries)
		cv_data.ids_obs = cv_train_entries
		run_checks(cv_data, raw_data)
		
		# loop thru hyperparams
		for (n, hps) in enumerate(hyperparams)
			# train model
			cv_model, cv_losses = construct_train_model(hps, cv_data, raw_data, nb_epochs, α=α, record_loss=record_loss, use_adam=use_adam)
            if record_loss
                _viz_loss!(axs[i_fold, n], cv_losses)
            end
			
			# evaluate on cv-test data
			perf_metric = compute_perf_metrics(cv_model, raw_data, cv_test_entries)
			perf_metrics[i_fold, n] = perf_metric[type_of_perf_metric]
		end		
	end

    mean_perf_over_folds = mean(eachrow(perf_metrics))[:]
    @assert length(mean_perf_over_folds) == length(hyperparams)
    opt_hps_id = argmax(mean_perf_over_folds)
	opt_hyperparams = hyperparams[opt_hps_id]

	return perf_metrics, opt_hps_id, opt_hyperparams, fig_losses
end
