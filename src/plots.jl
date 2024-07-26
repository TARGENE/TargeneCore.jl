function TMLE.pvalue(x)
    return try
        TMLE.pvalue(significance_test(x))
    catch
        NaN
    end
end

TMLE.pvalue(x::TargetedEstimation.FailedEstimate) = missing

log10_uniform_quantiles(n) = -log10.(collect(LinRange(0., 1., n + 1))[2:end])

log10_beta_quantiles(n, alpha) = -log10.([quantile(Beta(k, n + 1 − k), alpha) for k in 1:n])

function load_results(resultsfile; verbosity=0)
    verbosity > 0 && @info "Loading results."
    results = jldopen(io -> io["results"], resultsfile)
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

function summary_plots(results_file;
    outprefix="final",
    verbosity=0
    )
    results = load_results(results_file; verbosity=verbosity)
    qqplot(results, outprefix)
end