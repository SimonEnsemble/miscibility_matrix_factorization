function build_Xy(data::MiscibilityData, raw_data::RawData)
    X_train = []
    y_train = []
    X_test = []
    y_test = []
    # loop over ALL pairs. be sure to repeat. 
    # this is important to capture permutation-invariance
    for i = 1:data.n_compounds
        for j = 1:data.n_compounds
            if i == j
                continue
            end
            # build feature vector for this compound pair.
            x_ij = vcat(raw_data.X[:, i], raw_data.X[:, j])
            y_ij = raw_data.M_complete[i, j]
            
            if ! ismissing(data.M[i, j])
                push!(X_train, x_ij)
                push!(y_train, y_ij)
            else
                push!(X_test, x_ij)
                push!(y_test, y_ij)
            end
        end
    end
    @assert length(data.ids_obs) * 2 == length(y_train)
    @assert length(data.ids_missing) * 2 == length(y_test)
    @assert length(X_train) == length(y_train)
    @assert length(X_test) == length(y_test)
    @assert length(X_train[1]) == (7+4)*2 # number of features
    return X_train, y_train, X_test, y_test
end

function rf_feature_importance(
    data::MiscibilityData,
    raw_data::RawData
)
    X_train, y_train, X_test, y_test = build_Xy(data, raw_data)
    
    # train
    rf = RandomForestClassifier()
    rf.fit(X_train, y_train)

    res = permutation_importance(
        rf, X_test, y_test, n_repeats=100, scoring="balanced_accuracy"
    )
    
    # features replicated twice to handle permutation invariance
    nb_features = length(raw_data.features)
    mean_importance = (res["importances_mean"][1:nb_features] + res["importances_mean"][nb_features+1:end]) / 2
    std_importance = sqrt.((res["importances_std"][1:nb_features] .^ 2 + res["importances_std"][nb_features+1:end] .^ 2) / 2)
    return mean_importance, std_importance
end

function test_perf_baseline_model(
    data::MiscibilityData,
    raw_data::RawData;
    set_opt_cutoff::Bool=true
)
    X_train, y_train, X_test, y_test = build_Xy(data, raw_data)
    
    # train
    rf = RandomForestClassifier()
    rf.fit(X_train, y_train)

    # find opt cutoff
    if set_opt_cutoff
        prob_preds = rf.predict_proba(X_train)
        cutoffs = range(0.1, 0.9, length=20)
        scores = zeros(20)
        for (n, cutoff) in enumerate(cutoffs) 
            ŷ_train = [p > cutoff for p in prob_preds[:, 2]]
            scores[n] = balanced_accuracy_score(y_train, ŷ_train)
        end

        cutoff = cutoffs[argmax(scores)]
    else
        cutoff = 0.5
    end

    # test
    prob_preds = rf.predict_proba(X_test)
    ŷ_test = [p > cutoff for p in prob_preds[:, 2]]
    #ŷ_test = rf.predict(X_test)

    return (f1=f1_score(y_test, ŷ_test),
            accuracy=accuracy_score(y_test, ŷ_test),
            precision=precision_score(y_test, ŷ_test),
            recall=recall_score(y_test, ŷ_test),
            balanced_accuracy=balanced_accuracy_score(y_test, ŷ_test),
            cm=confusion_matrix(y_test, ŷ_test),
            cutoff=cutoff
    )
end

function miscibility_probabilities(data::MiscibilityData, raw_data::RawData)
    compound_categories = unique(raw_data.classes)

    π_miscible = Dict{String, Dict{String, Float64}}()
    for cat_a in compound_categories
        π_miscible[cat_a] = Dict{String, Float64}()
        for cat_b in compound_categories
            n_total = 0 # total # mixtures in this category
            n_miscible = 0
            for i = 1:raw_data.n_compounds, j = 1:raw_data.n_compounds
                if i == j
                    continue
                end
                if (cat_a, cat_b) == (raw_data.classes[i], raw_data.classes[j])
                    if ! ismissing(data.M[i, j])
                        n_total += 1
                        n_miscible += data.M[i, j]
                    end
                end
            end
            π_miscible[cat_a][cat_b] = n_miscible / n_total
        end
    end
    return π_miscible
end

function test_perf_guessing(data::MiscibilityData, raw_data::RawData)
    π_miscible = miscibility_probabilities(data, raw_data)
    
    m = [raw_data.M_complete[i, j] for (i, j) in data.ids_missing]
    m̂ = [NaN for i = 1:length(data.ids_missing)]
    for (k, (i, j)) in enumerate(data.ids_missing)
        m̂[k] = rand() < π_miscible[raw_data.classes[i]][raw_data.classes[j]]
    end

    return (f1=f1_score(m, m̂),
            accuracy=accuracy_score(m, m̂),
            precision=precision_score(m, m̂),
            recall=recall_score(m, m̂),
            balanced_accuracy=balanced_accuracy_score(m, m̂),
            cm=confusion_matrix(m, m̂)
    )
end
