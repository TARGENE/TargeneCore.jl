GRMIDs(file) = CSV.File(file, 
                        header=["FAMILY_ID", "SAMPLE_ID"], 
                        select=["SAMPLE_ID"],
                        types=String) |> DataFrame

function readGRM(prefix)
    ids = GRMIDs(string(prefix, ".id"))
    n = size(ids, 1)
    grm_size = n*(n + 1) ÷ 2
    GRM = mmap(string(prefix, ".bin"), Vector{Float32}, grm_size)

    return GRM, ids
end

function align_ic(ic, sample_ids, grm_ids)
    leftjoin!(
        grm_ids, 
        DataFrame(IC=ic, SAMPLE_ID=sample_ids),
        on=:SAMPLE_ID
    )
    aligned_ic = grm_ids.IC
    select!(grm_ids, Not(:IC))
    return coalesce.(aligned_ic, 0)
end

"""
    bit_distances(sample_grm, nτs)

Returns a matrix of shape (n_samples,nτs) where n_samples is the 
size of sample_grm.
The sample_grm comes from the gcta software. 
The process is as follows:
- Round between -1 and 1 as some elements may be beyond that limit
- Take 1 - this value to convert this quantity to a distance
- For each τ return if the distance between individuals is less or equal than τ
"""
function bit_distances(sample_grm, τs)
    distances = 1 .-  max.(min.(sample_grm, 1), -1)
    return convert(Matrix{Float32}, permutedims(distances) .<= τs)
end


default_τs(nτs;max_τ=2) = Float32[max_τ*(i-1)/(nτs-1) for i in 1:nτs]

retrieve_sample_ids(sample_ids::AbstractVector, batch_results) = sample_ids

retrieve_sample_ids(index::Int, batch_results) = batch_results[index].SAMPLE_IDS

function update_work_lists_with!(result::TMLE.JointEstimate, sample_ids, batch_results, grm_ids, results, influence_curves, n_obs)
    for estimate in result.estimates
        update_work_lists_with!(estimate, sample_ids, batch_results, grm_ids, results, influence_curves, n_obs)
    end
end

function update_work_lists_with!(result, sample_ids, batch_results, grm_ids, results, influence_curves, n_obs)
    if length(result.IC) > 0
        sample_ids = string.(retrieve_sample_ids(sample_ids, batch_results))
        push!(influence_curves, align_ic(result.IC, sample_ids, grm_ids))
        push!(n_obs, size(sample_ids, 1))
        push!(results, result)
    end
end

function build_work_list(prefix, grm_ids; estimator_key=1)
    dirname_, prefix_ = splitdir(prefix)
    dirname__ = dirname_ == "" ? "." : dirname_
    hdf5files = filter(
            x -> startswith(x, prefix_) && endswith(x, ".hdf5"), 
            readdir(dirname__)
    )
    hdf5files = sort([joinpath(dirname_, x) for x in hdf5files])

    influence_curves = Vector{Float32}[]
    n_obs = Int[]
    results = []
    for hdf5file in hdf5files
        jldopen(hdf5file) do io
            for key in keys(io)
                batch_results = io[key]
                for nt_result in batch_results
                    result = nt_result[estimator_key]
                    result isa TMLECLI.FailedEstimate && continue
                    sample_ids = nt_result.SAMPLE_IDS
                    update_work_lists_with!(
                        result,
                        sample_ids,
                        batch_results, 
                        grm_ids, results, 
                        influence_curves, 
                        n_obs
                    )
                end
            end
        end
    end
    influence_curves = length(influence_curves) > 0 ? reduce(vcat, transpose(influence_curves)) : Matrix{Float32}(undef, 0, 0)
    return results, influence_curves, n_obs
end


"""
    normalize(variances, n_observations)

Divides the variance estimates by the effective number of observations 
used for each phenotype at estimation time.
"""
normalize!(variances, n_observations) = 
    variances ./= permutedims(n_observations)


"""
    aggregate_variances(influence_curves, indicator, sample)

This function computes the sum for a single index i, see also `compute_variances`.
As the GRM is symetric it is performed as : 
    2 times off-diagonal elements with j < i + diagonal term 
and this for all τs. Some diagonal elements are not equal to 1 while in theory they should be,
this is an artefact of the GRM. We thus always include diagonal elements in the aggregate 
whatever the value of the indicator.
"""
function aggregate_variances(influence_curves, indicator, sample)
    @views begin
        D_off_diag = transpose(influence_curves[:, 1:sample-1])
        D_diag = transpose(influence_curves[:, sample])
        return D_diag .* (2indicator[:, 1:sample-1] * D_off_diag .+ D_diag)
    end
end


"""
    compute_variances(influence_curves, nτs, grm_files)

An overall variance estimate for a distance function d, a threshold τ 
and influence curve D is given by:
            σ̂ = 1/n ∑ᵢ∑ⱼ1(dᵢⱼ ≤ τ)DᵢDⱼ
              = 1/n 2∑ᵢDᵢ∑ⱼ<ᵢ1(dᵢⱼ ≤ τ)DᵢDⱼ + 1(dᵢᵢ ≤ τ)DᵢDᵢ

This function computes those variance estimates at each τ for all phenotypes
and queries.

# Arguments:
- influence_curves: Array of size (n_samples, n_queries, n_phenotypes)
- grm: Vector containing the lower elements of the GRM ie: n(n+1)/2 elements
- τs: 1 row matrix containing the distance thresholds between individuals
- n_obs: Vector of size (n_phenotypes,), containing the number of effective 
observations used during estimation

# Returns:
- variances: An Array of size (nτs, n_curves) where n_curves is the number of influence curves.
"""
function compute_variances(influence_curves, grm, τs, n_obs)
    n_curves, n_samples = size(influence_curves)
    variances = zeros(Float32, length(τs), n_curves)
    start_idx = 1
    for sample in 1:n_samples
        # lower diagonal of the GRM are stored in a single vector 
        # that are accessed one row at a time
        end_idx = start_idx + sample - 1
        sample_grm = view(grm, start_idx:end_idx)
        indicator = bit_distances(sample_grm, τs)
        variances .+= aggregate_variances(influence_curves, indicator, sample)
        start_idx = end_idx + 1
    end
    normalize!(variances, n_obs)
    return variances
end

function grm_rows_bounds(n_samples)
    bounds = Pair{Int, Int}[]
    start_idx = 1
    for sample in 1:n_samples
        # lower diagonal of the GRM are stored in a single vector 
        # that are accessed one row at a time
        end_idx = start_idx + sample - 1
        push!(bounds, start_idx => end_idx)
        start_idx = end_idx + 1
    end
    return bounds
end

function save_results(filename, results, τs, variances)
    jldopen(filename, "w") do io
        io["taus"] = τs
        io["variances"] = variances
        io["results"] = results
    end
end

corrected_stderrors(variances) =
    sqrt.(view(maximum(variances, dims=1), 1, :))

with_updated_std(estimate::T, std) where T <: TMLE.Estimate = T(
    estimate.estimand,
    estimate.estimate,
    convert(Float64, std),
    estimate.n,
    Float64[]
)

with_updated_std(results::AbstractVector, stds) =
    [NamedTuple{(:SVP,)}([with_updated_std(result, std)]) for (result, std) in zip(results, stds)]


"""
    sieve_variance_plateau(input_prefix;
        out="svp.hdf5",
        grm_prefix="GRM",
        verbosity=0, 
        n_estimators=10, 
        max_tau=0.8,
        estimator_key=1
    )

Sieve Variance Plateau CLI.

# Args

- `input-prefix`: Input prefix to HDF5 files generated by the tmle CLI.

# Options

- `-o, --out`: Output filename.
- `-g, --grm-prefix`: Prefix to the aggregated GRM.
- `-v, --verbosity`: Verbosity level.
- `-n, --n_estimators`: Number of variance estimators to build for each estimate. 
- `-m, --max_tau`: Maximum distance between any two individuals.
- `-e, --estimator-key`: Estimator to use to proceed with sieve variance correction.
"""
function sieve_variance_plateau(input_prefix::String;
    out::String="svp.hdf5",
    grm_prefix::String="GRM",
    verbosity::Int=0, 
    n_estimators::Int=10, 
    max_tau::Float64=0.8,
    estimator_key::String="1"
    )
    estimator_key_ = tryparse(Int, estimator_key)
    estimator_key = estimator_key_ === nothing ? Symbol(estimator_key) : estimator_key_
    τs = default_τs(n_estimators;max_τ=max_tau)
    grm, grm_ids = readGRM(grm_prefix)
    verbosity > 0 && @info "Preparing work list."
    results, influence_curves, n_obs = build_work_list(input_prefix, grm_ids, estimator_key=estimator_key)

    if length(influence_curves) > 0
        verbosity > 0 && @info "Computing variance estimates."
        variances = compute_variances(influence_curves, grm, τs, n_obs)
        std_errors = corrected_stderrors(variances)
        results = with_updated_std(results, std_errors)
    else
        variances = Float32[]
    end
    save_results(out, results, τs, variances)

    verbosity > 0 && @info "Done."
    return 0
end
