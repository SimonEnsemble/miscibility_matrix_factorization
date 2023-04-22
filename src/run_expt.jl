function run_experiment(data::MiscibilityData,
                        raw_data::RawData,
                        graph_regularization::Bool,
                        nb_hyperparams::Int=10;
                        nfolds::Int=3,
                        nb_epochs::Int=250
)
    # create list of hyperparams to explore.
	hyperparams_cv = [(
			k=rand([2, 3]), 
			γ=graph_regularization ? rand(Uniform(0, 0.1)) : 0.0,
			λ=rand(),
			σ=nothing, 
			use_features=false)
			   for _ = 1:nb_hyperparams]

    # conduct hyper-param optimization viz K-folds cross validation on training data
    perf_metrics, opt_hyperparams, fig_losses = do_hyperparam_optimization(
		data, hyperparams_cv, raw_data, nb_epochs=nb_epochs)

    # train deployment model on all training data with opt hyper-params
	model, losses = construct_train_model(opt_hyperparams, data, raw_data, nb_epochs,
		record_loss=false)

    set_opt_cutoff!(model, raw_data, data.ids_obs)

    # compute performance metrics on train and test
    perf_metrics = (test =compute_perf_metrics(model, raw_data, data.ids_missing),
                    train=compute_perf_metrics(model, raw_data, data.ids_obs)
	)

    return (data=data, model=model, perf_metrics=perf_metrics,
            opt_hyperparams=opt_hyperparams, fig_losses=fig_losses
    )
end

function run_experiments(θ::Float64, raw_data::RawData;
				  		 nruns::Int=3,
				         nb_hyperparams::Int=10
)

	perf_data = DataFrame()
	all_metrics = [:f1, :accuracy, :precision, :recall, :balanced_accuracy]

	for r = 1:nruns
		println("run #$r/$nruns")
		### simulate data collection
		# introduce missing values into the matrix.
		# important to keep data the same for comparison.
		data = sim_data_collection(θ, raw_data, weigh_classes=false)

		### train the models.
		graph_regularization = true
		res = run_experiment(data, raw_data, graph_regularization, nb_hyperparams)
		for met in all_metrics # this way for Gadfly to work
			push!(perf_data, (run=r,
                              model="GR-MF",
							  # this metric
				 			  score=res.perf_metrics.test[met],
						      metric=String(met),
							  # opt hyperparams
						      k=res.opt_hyperparams.k,
						      γ=res.opt_hyperparams.γ,
				 			  λ=res.opt_hyperparams.λ,
							  cutoff=res.model.cutoff)
			)
		end

		# turn off graph regularization
		graph_regularization = false
		res = run_experiment(data, raw_data, graph_regularization, nb_hyperparams)
		for met in all_metrics # this way for Gadfly to work
			push!(perf_data, (run=r,
                              model="MF",
							  # this metric
				 			  score=res.perf_metrics.test[met],
						      metric=String(met),
							  # opt hyperparams
						      k=res.opt_hyperparams.k,
						      γ=res.opt_hyperparams.γ,
				 			  λ=res.opt_hyperparams.λ,
							  cutoff=res.model.cutoff)
			)
		end

		# random forest
		baseline = test_perf_baseline_model(data, raw_data, set_opt_cutoff=true)
		for met in all_metrics # this way for Gadfly to work
			push!(perf_data, (run=r,
                              model="RF",
							  # this metric
				 			  score=baseline[met],
						      metric=String(met),
							  # opt hyperparams
						      k=NaN,
						      γ=NaN,
				 			  λ=NaN,
							  cutoff=baseline[:cutoff]), promote=true
			)
		end
	end
	return perf_data
end
