class_to_marker = Dict("Polymer"    => :circle,
                       "Protein"    => :rect,
                       "Surfactant" => :diamond,
                       "Salt"       => :cross)
class_to_color = Dict(zip(["Polymer", "Protein", "Salt", "Surfactant"],
                          [RGB(0.0, 0.74736935, 1.0),
                           RGB(0.83092886, 0.7934697, 0.22566332),
                           RGB(1.0, 0.426407, 0.681544),
                           RGB(0.0, 0.71687746, 0.55441886)])) # from gadfly
# a = RGB.(Gadfly.Scale.color_discrete_hue().f(4))
miscibility_colormap = reverse(ColorSchemes.:batlow)
label_to_color = Dict(0 => miscibility_colormap[1], 1 => miscibility_colormap[end])
label_to_string = Dict(0 => "immiscible", 1 => "miscible")

function viz_imputations(model::MFModel, data::MiscibilityData, raw_data::RawData)
    m     = [raw_data.M_complete[i, j] for (i, j) in data.ids_missing]
    m̂_raw = [pred_mᵢⱼ(model, i, j)     for (i, j) in data.ids_missing]

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="σ(cᵢ ⋅ cⱼ + b)", ylabel="# pairs of solutions")
    for c in [0, 1]
        hist!(m̂_raw[m .== c], color=(label_to_color[c], 0.5), label=label_to_string[c], 
              strokewidth=1, linewidth=2, strokecolor="black", bins=range(0.0, 1.0, length=11))
    end
    vlines!(ax, model.cutoff, linestyle=:dash, color="black", linewidth=2)
    axislegend(position=:lt)
    ylims!(0, nothing)
    xlims!(0, 1)
    save("viz_imputations.pdf", fig)
    return fig
end

function pretty_compound_labels(raw_data::RawData)
    the_compound_labels = rich.(raw_data.compounds)
    # subscripts
    the_compound_labels[findfirst(raw_data.compounds .== "K2HPO4")] = rich("K", subscript("2"), "HPO", subscript("4"))
    the_compound_labels[findfirst(raw_data.compounds .== "Na2HPO4")] = rich("Na", subscript("2"), "HPO", subscript("4"))
    return the_compound_labels
end

function viz_miscibility_matrix(M, raw_data::RawData; 
        draw_brackets::Bool=false, 
        include_legend::Bool=true, 
        savename::Union{String, Nothing}=nothing,
        show_solute_labels::Bool=true,
        show_xy_labels::Bool=true,
        huge_font::Bool=false
)
    big_fontsize = 50
    if huge_font
        big_fontsize *= 1.75
    end
    
    the_compound_labels = pretty_compound_labels(raw_data)

    fig = Figure(resolution=(1450, 1450))
    ax  = Axis(fig[1, 1], 
        xlabel=show_xy_labels ? "solution" : "", ylabel=show_xy_labels ? "solution" : "", 
        xgridvisible=false, 
        xlabelsize=big_fontsize,
        ylabelsize=big_fontsize,
        ygridvisible=false, 
        xticks=(1:raw_data.n_compounds, show_solute_labels ? the_compound_labels : ["" for _ = 1:raw_data.n_compounds]),
        yticks=(1:raw_data.n_compounds, show_solute_labels ? reverse(the_compound_labels) : ["" for _ = 1:raw_data.n_compounds]),
        xticklabelrotation=π/2
    )
    if ! draw_brackets
        ax.aspect = DataAspect()
    end
    heatmap!(ax, reverse(M', dims=2), colormap=miscibility_colormap)
    for i = 1:raw_data.n_compounds+1
        hlines!(ax, i - 0.5, color="gray", linewidth=1)
        vlines!(ax, i - 0.5, color="gray", linewidth=1)
    end
    
    if draw_brackets
        # add brackets to indicate the class of the compound
        t_ax = Axis(fig[1, 1, Top()], height=50, xautolimitmargin=(0, 0))
        r_ax = Axis(fig[1, 1, Right()], width=50, yautolimitmargin=(0, 0))
        for a in [t_ax, r_ax]
            hidedecorations!(a)
            hidespines!(a)
        end

        c0 = 0.5
        for c in unique(raw_data.classes) # loop thru classes
            l = sum(raw_data.classes .== c) # number of instances of this class

            # draw brackets on top
            bracket!(t_ax, c0 + 0.5, 0, c0 + l - 0.5, 0, orientation=:up, fontsize=big_fontsize/1.8,
                     linewidth=huge_font ? 8 : 5,
                     font=AlgebraOfGraphics.firasans("Light"), text=huge_font ? "" : lowercase(c), color=class_to_color[c])

            # draw brackets on right
            bracket!(r_ax, 0, raw_data.n_compounds + 1 - (c0 + l) + 0.5, 0, raw_data.n_compounds + 1 - c0 - 0.5, 
                     orientation=:down, fontsize=big_fontsize/1.8,
                     linewidth=huge_font ? 8 : 5,
                     font=AlgebraOfGraphics.firasans("Light"), text=huge_font ? "" : lowercase(c), color=class_to_color[c])
            
            c0 += l
        end

        linkxaxes!(ax, t_ax)
        linkyaxes!(ax, r_ax)
        ylims!(t_ax, 0, 2) # to see text
        xlims!(r_ax, 0, 2) # to see text
        colsize!(fig.layout, 1, Aspect(1, 1.0))
    end
    resize_to_layout!(fig)

    if include_legend
        legend_patches = [
            PolyElement(color=miscibility_colormap[1], strokecolor="gray", polystrokewidth=1), 
            PolyElement(color=miscibility_colormap[end], strokecolor="gray", polystrokewidth=1)
        ]
        legend_labels = ["immiscible", "miscible"]
        legend_patchsize = (35, 35)
        if any(ismissing.(M))
            push!(legend_patches, PolyElement(color="white", strokecolor="gray", polystrokewidth=1))
            push!(legend_labels, "missing")
            legend_patchsize = (35, 35, 35)
        end
        Legend(fig[0, 1], legend_patches, legend_labels, patchsize=legend_patchsize, orientation=:horizontal, labelsize=big_fontsize)
    end

    if draw_brackets
        notify(t_ax.finallimits)
        notify(r_ax.finallimits)
    end
    if ! isnothing(savename)
        if split(savename, ".")[end] != "pdf"
            savename *= ".pdf"
        end
        save(savename, fig)
    end
    return fig
end

function viz_C(
        model::MFModel, 
        raw_data::RawData; 
        draw_brackets::Bool=true, 
        savename::Union{String, Nothing}=nothing,
        minimal_viz::Bool=false
)
    big_fontsize = 50

    the_compound_labels = pretty_compound_labels(raw_data)
    colormap = reverse(ColorSchemes.diverging_gwr_55_95_c38_n256)

    fig = Figure(resolution=(1450, 1450))#(size(model.C)[1] * 1450 / raw_data.n_compounds, 1450))
    ax  = Axis(fig[1, 1],
        xlabel="solution",
        ylabel="latent\ndim.",
        xgridvisible=false,
        xlabelsize=big_fontsize,
        ylabelsize=big_fontsize,
        ygridvisible=false,
        xticks=(1:raw_data.n_compounds, the_compound_labels),
        yticks=(1:size(model.C)[1], ["$i" for i in reverse(1:size(model.C)[1])]),
        xticklabelrotation=π/2,
        aspect=DataAspect()
    )
    max_value = maximum(abs.(model.C))
    hm = heatmap!(ax, reverse(model.C', dims=2), colormap=colormap, colorrange=(-max_value, max_value))
    for i = 1:raw_data.n_compounds+1
        vlines!(ax, i - 0.5, color="gray", linewidth=1)
    end
    for i = 1:size(model.C)[1]+1
        hlines!(ax, i - 0.5, color="gray", linewidth=1)
    end
    if draw_brackets
        # add brackets to indicate the class of the compound
        t_ax = Axis(fig[1, 1, Top()], height=50, xautolimitmargin=(0, 0))
        hidedecorations!(t_ax)
        hidespines!(t_ax)

        c0 = 0.5
        for c in unique(raw_data.classes) # loop thru classes
            l = sum(raw_data.classes .== c) # number of instances of this class

            # draw brackets on top
            bracket!(t_ax, c0 + 0.5, 0, c0 + l - 0.5, 0, orientation=:up, fontsize=big_fontsize/1.8, 
                     linewidth=minimal_viz ? 8 : 5,
                     font=AlgebraOfGraphics.firasans("Light"), text=minimal_viz ? "" : lowercase(c), color=class_to_color[c])

            c0 += l
        end

        linkxaxes!(ax, t_ax)
        ylims!(t_ax, 0, 2) # to see text
    end
    #colsize!(fig.layout, 1, Aspect(1, 1.0))
    rowsize!(fig.layout, 1, Fixed(pixelarea(ax.scene)[].widths[2]))
    if ! minimal_viz
        Colorbar(fig[:, end+1], hm, label="Cᵢⱼ")
    end
    # this messes up the brackets
    if minimal_viz
        hidedecorations!(ax)
    end
    resize_to_layout!(fig)
    if draw_brackets
       notify(t_ax.finallimits)
       notify(ax.finallimits)
    end
    if ! isnothing(savename)
        save(savename, fig)
    end
    return fig
end

_viz_loss!(ax, losses::Vector{Float64}) = lines!(ax, 1:length(losses), losses)

function viz_loss(losses::Vector{Float64}; save_fig::Bool=false, append_filename::String="")
    fig = Figure(backgroundcolor=("white", 1.0))
    ax = Axis(fig[1, 1], xlabel="# epochs", ylabel="loss")
    _viz_loss!(ax, losses)
    if save_fig
        save("loss" * append_filename * ".pdf", fig)
    end
    fig
end

function viz_latent_space_3d(model::MFModel, raw_data::RawData; append_filename::String="")
    if size(model.C)[1] != 3
        error("this for k = 3")
    end

    fig = Figure()
    ax = Axis3(fig[1, 1],
        xlabel="latent dim. 1",
        ylabel="latent dim. 2",
        zlabel="latent dim. 3",
        xgridvisible=true,
        ygridvisible=true,
        zgridvisible=true,
        aspect=:data
    )
    # draw axes
    lines!([-5, 5], [0, 0], [0, 0], color="black", linewidth=1)
    lines!([0, 0], [-5, 5], [0, 0], color="black", linewidth=1)
    lines!([0, 0], [0, 0], [-5, 5], color="black", linewidth=1)
    classes = unique(raw_data.classes)
    sps = []
    for c in classes
        push!(sps,
            scatter!(model.C[:, raw_data.classes .== c],
                color=class_to_color[c], label=c)
        )
    end
    the_lims = 1.1 * maximum(abs.(model.C))
    xlims!(-the_lims, the_lims)
    ylims!(-the_lims, the_lims)
    zlims!(-the_lims, the_lims)
    Legend(fig[1, 2], sps, lowercase.(classes), orientation=:vertical)
    #resize_to_layout!(fig)
    save("3d_latent_space" * append_filename * ".pdf", fig)
    return fig
end

function viz_latent_space(model::MFModel, raw_data::RawData; incl_legend::Bool=true, save_fig::Bool=false, append_filename::String="")
    do_pca = size(model.C)[1] != 2

    if do_pca
        # standardize for PCA
        C̃ = deepcopy(model.C')
        for c = 1:size(C̃)[2]
            C̃[:, c] = (C̃[:, c] .- mean(C̃[:, c])) / std(C̃[:, c])
        end
        # the projection into 2D via PCA
        pca = PCA(n_components=2)
        Ĉ = pca.fit_transform(C̃)'
    else
        Ĉ = model.C
    end

    fig = Figure()
    ax  = Axis(fig[1, 1],
        xlabel=do_pca ? "PC1" : "latent dimension 1",
        ylabel=do_pca ? "PC2" : "latent dimension 2",
        aspect=DataAspect()
    )

    vlines!(0.0, color="lightgray", linewidth=1, linestyle=:dash)
    hlines!(0.0, color="lightgray", linewidth=1, linestyle=:dash)
    sps = []
    classes = ["Polymer", "Protein", "Surfactant", "Salt"]
    for class in classes
        ids = raw_data.classes .== class
        push!(sps,
              scatter!(Ĉ[1, ids], Ĉ[2, ids],
                strokewidth=1, strokecolor="black", marker=class_to_marker[class],
                color=class_to_color[class], label=class, aspect=DataAspect())
             )
    end
    if incl_legend
        Legend(fig[0, 1], sps, lowercase.(classes), orientation=:horizontal)
    end
    if save_fig
        save("latent_space" * append_filename * ".pdf", fig)
    end
    return fig
end

function viz_confusion(cm::Matrix; save_fig::Bool=false)
    cm_to_plot = reverse(cm, dims=1)'

    fig = Figure(resolution=(300, 300))
    ax  = Axis(fig[1, 1],
        xlabel="prediction", ylabel="truth",
        xticks=(1:2, ["immiscible", "miscible"]),
        aspect=DataAspect(),
        yticks=(1:2, reverse(["immiscible", "miscible"])),
        yticklabelrotation=π/2
    )

    # make negative to use negative region of colormap
    cm_to_plot[1, 2] = -cm_to_plot[1, 2]
    cm_to_plot[2, 1] = -cm_to_plot[2, 1]
    hm = heatmap!(cm_to_plot,
        colormap=ColorSchemes.diverging_gwr_55_95_c38_n256,
        colorrange=(-maximum(abs.(cm)), maximum(abs.(cm)))
    )
    
    cm_to_plot[1, 2] = -cm_to_plot[1, 2]
    cm_to_plot[2, 1] = -cm_to_plot[2, 1]
    for i = 1:2
        for j = 1:2
            text!("$(round(Int, cm_to_plot[i, j]))",
                  position=(i, j), align=(:center, :center), color="black",
                  fontsize=45
            )
        end
    end
    #Colorbar(fig[1, 2], hm, label="# pairs")
    resize_to_layout!(fig)
    if save_fig
        save("confusion_matrix.pdf", fig)
    end
    return fig
end

function viz_category_miscibility(raw_data::RawData)
    compound_categories = unique(raw_data.classes)

    M_categories = zeros(length(compound_categories), length(compound_categories))

    Ω_all = [(i, j) for i = 1:raw_data.n_compounds for j = i+1:raw_data.n_compounds]
    @assert length(Ω_all) == Int(raw_data.n_compounds*(raw_data.n_compounds-1)/2)

    for c = 1:4
        for c′ = 1:4
            function pair_in_category_pair(i, j)
                this_cat = (raw_data.classes[i], raw_data.classes[j])

                return ((this_cat == (compound_categories[c], compound_categories[c′])) || (this_cat == (compound_categories[c′], compound_categories[c])))
            end
            # get pairs belonging to these two categories
            ids = [(i, j) for (i, j) in Ω_all if pair_in_category_pair(i, j)]

            outcomes = [raw_data.M_complete[i, j] for (i, j) in ids]

            M_categories[c, c′] = 1.0 - mean(outcomes)
            if c > c′
                M_categories[c, c′] = NaN
            end
        end
    end

    M_categories_plt = reverse(M_categories, dims=1)'

    fig = Figure()
    ax  = Axis(fig[1, 1],
        xlabel="category", ylabel="category",
        xticks=(1:4, lowercase.(compound_categories)),
        aspect=DataAspect(),
        yticks=(1:4, reverse(lowercase.(compound_categories))),
        yticklabelrotation=π/2
    )
    hm = heatmap!(M_categories_plt,
                  colormap=reverse(ColorSchemes.viridis), colorrange=(0, 1))

    for i = 1:4
        for j = 1:4
            if ! isnan(M_categories_plt[i, j])
                text!("$(round(M_categories_plt[i, j], digits=2))",
                      position=(i, j), align=(:center, :center), color="black",
                      fontsize=25
                )
            end
        end
    end
    Colorbar(fig[1, 2], hm, label="fraction ATPS forming")
    save("category_based_miscibility.pdf", fig)
    return fig
end

function viz_rf_feature_importance(
    raw_data::RawData,
    θ::Float64,
    nb_runs::Int
)
    # get feature importance
    μ_importance, σ_importance = rf_feature_importance(raw_data, θ, nb_runs)

    fig = Figure()
    ax  = Axis(fig[1, 1],
        xlabel="solution feature", ylabel="decrease of\nbalanced accuracy",
        xticks=(1:length(raw_data.features), raw_data.features),
        xticklabelrotation=π/2,
        title="permutation-based feature importance"
    )
    bar_colors = map(i -> i > 0 ? "green" : "red", μ_importance)
    barplot!(1:length(raw_data.features), μ_importance, color=bar_colors)
    errorbars!(1:length(raw_data.features), μ_importance,
        σ_importance, whiskerwidth=10)
    Label(fig[1, 1], "θ = $(round(θ, digits=1))", tellwidth=false, tellheight=false, valign=0.9, halign=0.9)
    save("rf_feature_importance.pdf", fig)
    fig
end

function viz_ml_procedure(data::MiscibilityData, raw_data::RawData, model::MFModel; nfolds::Int=3)
    if ! isdir("viz_procedure")
        mkdir("viz_procedure")
    end
    
    viz_C(model, raw_data, savename="viz_procedure/C.pdf", draw_brackets=false, minimal_viz=true)
    hack_model = deepcopy(model)
    fill!(hack_model.C, hack_model.b)
    viz_C(hack_model, raw_data, savename="viz_procedure/b.pdf", draw_brackets=false, minimal_viz=true)

    # viz complete matrix
    viz_miscibility_matrix(raw_data.M_complete, raw_data, 
        draw_brackets=true, show_solute_labels=false,
        savename="viz_procedure/M_complete.pdf", include_legend=true, huge_font=true)

    viz_kwargs = (draw_brackets=false, show_solute_labels=false, show_xy_labels=false, include_legend=false, huge_font=true)
    # viz training data
    viz_miscibility_matrix(data.M, raw_data; viz_kwargs..., savename="viz_procedure/M_train.pdf")
    viz_miscibility_matrix(data.M, raw_data; viz_kwargs..., savename="viz_procedure/M_train_w_label.pdf", show_xy_labels=true, include_legend=true)

    # build test data
    M_test = deepcopy(data.M)
    fill!(M_test, missing)
    for (i, j) in data.ids_missing
        M_test[i, j] = M_test[j, i] = raw_data.M_complete[i, j]
    end
    viz_miscibility_matrix(M_test, raw_data; viz_kwargs..., savename="viz_procedure/M_test.pdf", show_xy_labels=true)

    # 3-folds cross-validation split
    kf_split = train_test_pairs(StratifiedCV(nfolds=nfolds, shuffle=true),
                            1:length(data.ids_obs),
                            [data.M[i, j] for (i, j) in data.ids_obs])

    for (id_fold, (ids_cv_train, ids_cv_test)) in enumerate(kf_split)
        # get the list of matrix entries, vector of tuples
        cv_train_entries = data.ids_obs[ids_cv_train]
        cv_test_entries  = data.ids_obs[ids_cv_test]

        ###
        # create copy of data
        # introduce additional missing values, where the cv-test data are.
        M_train = deepcopy(data.M)
        for (i, j) in cv_test_entries
            # ablate entries
            M_train[i, j] = missing
            M_train[j, i] = missing
        end
        M_test = deepcopy(data.M)
        for (i, j) in cv_train_entries
            M_test[i, j] = missing
            M_test[j, i] = missing
        end
        for i = 1:raw_data.n_compounds
            M_test[i, i] = missing
        end
        viz_miscibility_matrix(M_test, raw_data; viz_kwargs...,
            savename="viz_procedure/M_fold_$(id_fold)_cv_test.pdf"
           )
        viz_miscibility_matrix(M_train, raw_data;  viz_kwargs...,
            savename="viz_procedure/M_fold_$(id_fold)_cv_train.pdf"
           )
    end
    # the imputed matrix, too.
    M_predicted = pred_M(model) .> model.cutoff
    viz_miscibility_matrix(M_predicted, raw_data; savename="viz_procedure/M_imputed.pdf", viz_kwargs...)
    return "see viz_procedure/ for images"
end

function viz_hyperparams(hps, opt_hps)
    fig = Figure()
    ax = Axis3(fig[1, 1], 
               xlabel="k", ylabel="λ", zlabel="γ", 
               xticks=[2, 3], yticks=[0.2, 1.0], zticks=[0, 0.1],
               title="hyperparameter space")
    scatter!(
            [hp[:k] for hp in hps],
            [hp[:λ] for hp in hps],
            [hp[:γ] for hp in hps]
    )
    scatter!([opt_hps[:k]], [opt_hps[:λ]], [opt_hps[:γ]], markersize=15)
    save("hp_space.png", fig,  px_per_unit=2)
    return fig
end
