function add_pvalue_columns!(results)
    for colname ∈ names(results)
        results[!, Symbol(colname, :_PVALUE)] = [pvalue_or_nan(Ψ̂) for Ψ̂ ∈ results[!, colname]]
    end
end

TMLE.estimate(Ψ̂::TMLECLI.FailedEstimate) = "Failed"

get_treatment_changes(Ψ) = Dict(treatment_name => join(treatment_values, " => ") for (treatment_name, treatment_values) ∈ Ψ.treatment_values)
get_treatment_changes(Ψ::TMLE.JointEstimand) = [get_treatment_changes(Ψᵢ) for Ψᵢ ∈ Ψ.args]

get_estimand_type(Ψ) = replace(string(typeof(Ψ)), "TMLE.Statistical" => "")
get_estimand_type(Ψ::TMLE.JointEstimand) = get_estimand_type(first(Ψ.args))

log10_uniform_quantiles(n) = -log10.(collect(LinRange(0., 1., n + 1))[2:end])

log10_beta_quantiles(n, alpha) = -log10.([quantile(Beta(k, n + 1 − k), alpha) for k in 1:n])

function qqplot(output, results)
    fig = Figure()
    ax = Axis(fig[1, 1], 
        title="Q-Q Plot",
        xlabel = "-log₁₀(Uniform)",
        ylabel = "-log₁₀(Observed)"
    )
    # QQ plots
    pvalue_cols = [colname for colname ∈ names(results) if endswith(colname, "PVALUE")]
    for pvalue_col ∈ pvalue_cols
        pvalues = -log10.(filter(x -> !isnan(x), results[!, pvalue_col]))
        n = length(pvalues)
        unif_quantiles = log10_uniform_quantiles(n)
        qqplot!(ax, unif_quantiles, pvalues, qqline=:identity, label=replace(string(pvalue_col), "_PVALUE" => ""))
    end
    # Confidence bands (from: https://lbelzile.github.io/lineaRmodels/qqplot.html)
    n = size(results, 1)
    unif_quantiles = log10_uniform_quantiles(n)
    ub = log10_beta_quantiles(n, 0.025)
    lb = log10_beta_quantiles(n, 0.975)
    band!(ax, unif_quantiles, lb, ub, color=(:grey, 0.6))
    # Legend
    axislegend(position=:lt)
    # Save figure
    CairoMakie.save(output, fig)
    return fig
end

function save_summary_yaml(output, results)
    estimator_names = Symbol.(filter(x -> !endswith(x, "PVALUE"), names(results)))
    estimator_1 = first(estimator_names)
    summary_dicts = []
    for row in eachrow(results)
        Ψ = row[estimator_1].estimand
        estimand_results = OrderedDict(
            :EFFECT_TYPE => get_estimand_type(Ψ),
            :OUTCOME => get_outcome(Ψ),
            :TREATMENTS => get_treatment_changes(Ψ),            
        )
        for estimator_name in estimator_names
            estimand_results[estimator_name] = OrderedDict(
                :PVALUE => row[Symbol(estimator_name, :_PVALUE)],
                :EFFECT_SIZE => TMLE.estimate(row[estimator_name])
            )
        end
        push!(summary_dicts, estimand_results)
    end
    YAML.write_file(output, summary_dicts)
end

function load_results_as_df(input_prefix)
    input_files = files_matching_prefix(input_prefix)
    results = DataFrame(TMLECLI.read_results_from_files(input_files))
    add_pvalue_columns!(results)
    return results
end

function make_outputs(input_prefix;
    output_prefix="final",
    verbosity=1
    )
    # Load results
    verbosity > 0 && @info("Loading results from files.")
    results = load_results_as_df(input_prefix)
    # Write full results
    verbosity > 0 && @info("Saving aggregated results.")
    jldsave(make_filepath_from_prefix(output_prefix; filename="results.hdf5"); results)
    # Write Summary YAML file
    verbosity > 0 && @info("Saving summary YAML file.")
    save_summary_yaml(make_filepath_from_prefix(output_prefix; filename="results.summary.yaml"), results)
    # Save QQ plot
    verbosity > 0 && @info("Making QQ plot.")
    qqplot(make_filepath_from_prefix(output_prefix; filename="QQ.png"), results)
    
    verbosity > 0 && @info("Done.")
    return 0
end