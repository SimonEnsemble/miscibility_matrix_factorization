module MiscibilityMF

using CSV, DataFrames, Statistics, LogExpFunctions, LinearAlgebra, Random, 
      StatsBase, ScikitLearn, Distributions, KernelFunctions, CairoMakie, ColorSchemes

import MLJBase: partition, StratifiedCV, train_test_pairs

import AlgebraOfGraphics
AlgebraOfGraphics.set_aog_theme!(fonts=[AlgebraOfGraphics.firasans("Light"), AlgebraOfGraphics.firasans("Light")])
update_theme!(fontsize=16)

@sk_import metrics: (confusion_matrix, accuracy_score, f1_score, precision_score, recall_score)
@sk_import ensemble: RandomForestClassifier
@sk_import decomposition: PCA

include.(
    [
     "data_prep.jl", "viz.jl", "introduce_missing_entries.jl", "kernel_matrix.jl", "grmf.jl", "perf.jl", "hyperparam_opt.jl"
    ]
)

export build_miscibility_matrix, build_compound_info, retreive_raw_data, RawData, # data_prep.jl
       viz_miscibility_matrix, viz_loss, viz_latent_space, # viz.jl
       sim_data_collection, # introduce_missing_entries.jl
       class_similarity_matrix, kernel_matrix, # kernel_matrix.jl
       construct_train_model, MFModel, # grmf.jl
       compute_perf_metrics, # perf.jl
       do_hyperparam_optimization # hyperparam_opt.jl
end
