module MiscibilityMF

using CSV, DataFrames, Statistics, LogExpFunctions, LinearAlgebra, Random, 
      StatsBase, ScikitLearn, Distributions, KernelFunctions, CairoMakie, ColorSchemes

import MLJBase: partition, StratifiedCV, train_test_pairs

@sk_import metrics: (confusion_matrix, accuracy_score, f1_score, precision_score, recall_score)
@sk_import ensemble: RandomForestClassifier
@sk_import decomposition: PCA

include.(
    [
     "data_prep.jl", "viz.jl", "introduce_missing_entries.jl", "kernel_matrix.jl", "grmf.jl"
    ]
)

export build_miscibility_matrix, build_compound_info, retreive_raw_data, # data_prep.jl
       viz_miscibility_matrix, # viz.jl
       sim_data_collection, # introduce_missing_entries.jl
       class_similarity_matrix, kernel_matrix, # kernel_matrix.jl
    construct_train_model # grmf.jl

end
