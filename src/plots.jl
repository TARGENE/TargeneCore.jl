TMLE.pvalue(x) = pvalue(significance_test(x))
TMLE.pvalue(x::TargetedEstimation.FailedEstimate) = missing

log10_uniform_quantiles(n) = -log10.(collect(LinRange(0., 1., n + 1))[2:end])

log10_beta_quantiles(n, alpha) = -log10.([quantile(Beta(k, n + 1 − k), alpha) for k in 1:n])

function read_results_file(file)
    jldopen(file) do io
        return reduce(vcat, (io[key] for key in keys(io)))
    end
end

function read_results_files(prefix)
    directory_, prefix_ = splitdir(prefix)
    directory = directory_ == "" ? "." : directory_
    matching_files = [joinpath(directory_, file) for file in readdir(directory) if startswith(file, prefix_)]
    reduce(vcat, (read_results_file(file) for file in matching_files))
end

function load_results(file_or_prefix; verbosity=0)
    verbosity > 0 && @info "Loading results."
    results = isfile(file_or_prefix) ? read_results_file(file_or_prefix) : read_results_files(file_or_prefix)
    estimators = collect(key for key ∈ keys(first(results)) if key !== :SAMPLE_IDS)
    results_df = DataFrame([[r[id] for r in results] for id in 1:length(estimators)], estimators)
    for estimator in estimators
        results_df[!, Symbol(estimator, :_PVALUE)] = [pvalue(x) for x in results_df[!, estimator]]
    end
    return results_df
end

function qqplot(results, outprefix)
    fig = Figure()
    ax = Axis(fig[1, 1], 
        title="Q-Q Plot",
        xlabel = "-log₁₀(Uniform)",
        ylabel = "-log₁₀(Observed)"
    )
    # QQ plots
    pvalue_cols = [colname for colname ∈ names(results) if endswith(colname, "PVALUE")]
    for pvalue_col ∈ pvalue_cols
        pvalues = -log10.(filter(x -> !isnan(x), skipmissing(results[!, pvalue_col])))
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
    CairoMakie.save(TargeneCore.make_filepath_from_prefix(outprefix; filename="QQ.png"), fig)
    return fig
end

function summary_plots(results_prefix;
    outprefix="final",
    verbosity=0
    )
    results = load_results(results_prefix; verbosity=verbosity)
    qqplot(results, outprefix)
end