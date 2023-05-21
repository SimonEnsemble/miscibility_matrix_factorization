module MiscibilityMF

using CSV, DataFrames, Statistics, LogExpFunctions, LinearAlgebra, Random, 
      StatsBase, ScikitLearn, Distributions, KernelFunctions, CairoMakie, ColorSchemes, Colors

import MLJBase: partition, StratifiedCV, train_test_pairs

import AlgebraOfGraphics
AlgebraOfGraphics.set_aog_theme!(fonts=[AlgebraOfGraphics.firasans("Light"), AlgebraOfGraphics.firasans("Light")])
update_theme!(fontsize=16)

# https://discourse.julialang.org/t/error-when-calling-scikitlearn-from-a-function-in-a-package/46131/8
function __init__()
    @eval @sk_import metrics: (confusion_matrix, accuracy_score, f1_score, precision_score, recall_score, balanced_accuracy_score)
    @eval @sk_import ensemble: RandomForestClassifier
    @eval @sk_import decomposition: PCA
end

include.(
    [
     "data_prep.jl", "introduce_missing_entries.jl", "kernel_matrix.jl", "grmf.jl", "perf.jl", "hyperparam_opt.jl", "set_cutoff.jl", "viz.jl", "baseline.jl", "run_expt.jl"
    ]
)

export build_miscibility_matrix, build_compound_info, retreive_raw_data, RawData, # data_prep.jl
       sim_data_collection, MiscibilityData, # introduce_missing_entries.jl
       class_similarity_matrix, kernel_matrix, # kernel_matrix.jl
       construct_train_model, MFModel, # grmf.jl
       compute_perf_metrics, # perf.jl
       do_hyperparam_optimization, # hyperparam_opt.jl
       set_opt_cutoff!, # set_cutoff.jl
       viz_miscibility_matrix, viz_loss, viz_latent_space, viz_confusion, class_to_color, # viz.jl
       build_Xy, test_perf_baseline_model, test_perf_guessing, # baseline.jl
       gen_hyperparams, run_experiment, run_experiments # run_expt.jl
end
