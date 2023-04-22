function run_experiment(θ::Float64,
                        raw_data::RawData,
                        use_features::Bool,
                        nhyperparams::Int=10, 
                        nfolds::Int=3, 
					    nb_epochs::Int=1000,
                        graph_regularization::Bool=true,
                        class_kernel::Bool=true,
	                    type_of_perf_metric::Symbol=:f1,
					    seed::Union{Nothing, Int}=nothing
)
	# create list of hyperparams to explore.
	if ! isnothing(seed)
		Random.seed!(seed)
	end
    data = sim_data_collection(θ, raw_data)
	
	hyperparams = [
		(k = rand([2, 3]),
		 # turn off graph reg on/off
		 γ = graph_regularization ? rand(Uniform(0, 0.1)) : 0.0, 
		 λ = rand(), 
		 # kernel based on class vs. features.
		 σ = class_kernel ? 10.0 : rand(Uniform(1e-1, 2))
		)
		    for _ = 1:ngrid]

	# conduct hyper-param optimization viz K-folds cross validation on training data
	perf_metric, opt_hyperparams, fig_losses = do_hyperparam_optimization(
		          data, cv_hyperparams, nfolds=nfolds, ngrid=ngrid, 
					    nb_epochs=nb_epochs, type_of_perf_metric = type_of_perf_metric)
	

	# train deployment model on all training data with opt hyper-params
	opt_model, opt_losses = construct_train_model(opt_hyperparams, data, nb_epochs)

	set_opt_cutoff!(opt_model, data, data.ids_obs)

	# compute performance metrics on train and test
	perf_metrics = (test =compute_perf_metrics(opt_model, data, data.ids_missing),
		            train=compute_perf_metrics(opt_model, data, data.ids_obs))

	# confusion matrix on test
	cm = compute_cm(opt_model, data)

	return (data=data, opt_model=opt_model, perf_metrics=perf_metrics, 
		    opt_hyperparams=opt_hyperparams, fig_losses=fig_losses, cm=cm
	)
end
