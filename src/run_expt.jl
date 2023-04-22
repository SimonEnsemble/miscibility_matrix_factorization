function run_experiment(θ::Float64,
                        raw_data::RawData,
                        graph_regularization::Bool,
                        nb_hyperparams::Int=10;
                        nfolds::Int=3,
                        nb_epochs::Int=250,
                        type_of_perf_metric::Symbol=:f1,
                        seed::Union{Nothing, Int}=nothing
)
    # create list of hyperparams to explore.
    if ! isnothing(seed)
        Random.seed!(seed)
    end
    data = sim_data_collection(θ, raw_data)

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

    baseline = test_perf_baseline_model(data, raw_data)

    # compute performance metrics on train and test
    perf_metrics = (test =compute_perf_metrics(model, raw_data, data.ids_missing),
                    train=compute_perf_metrics(model, raw_data, data.ids_obs),
	                baseline=baseline
	)

    return (data=data, model=model, perf_metrics=perf_metrics,
            opt_hyperparams=opt_hyperparams, fig_losses=fig_losses
    )
end
