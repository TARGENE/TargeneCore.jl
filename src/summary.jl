init_sieve(sieve_filename_::Nothing) = DataFrame()
init_sieve(sieve_filename_) =  DataFrame(
    SIEVE_STDERR = Union{Float64, Missing}[],
    SIEVE_PVAL = Union{Float64, Missing}[],
    SIEVE_LWB = Union{Float64, Missing}[],
    SIEVE_UPB = Union{Float64, Missing}[]
)

init_results(sieve_filename_) = hcat(DataFrame(
    PHENOTYPE = String[],
    PHENOTYPE_TYPE = String[],
    QUERYNAME = String[],
    FILENAME_ORIGIN = String[],
    REPORT_KEY = String[],
    INITIAL_ESTIMATE = Float64[],
    ESTIMATE = Float64[],
    STDERR = Float64[],
    PVAL = Float64[],
    LWB = Float64[],
    UPB = Float64[],
    ),
    init_sieve(sieve_filename_)
)

sieve_results_(estimate, sieve_stderror::Nothing) =
    [missing, missing, missing, missing]

function sieve_results_(estimate, sieve_stderror)
    sieve_test = OneSampleZTest(estimate, sieve_stderror, 1)
    sieve_pval = pvalue(sieve_test)
    sieve_lwb, sieve_upb = confint(sieve_test)
    return [sieve_stderror, sieve_pval, sieve_lwb, sieve_upb]
end

sieve_results(estimate, filename, key, pair_to_var_id::Nothing, sieve_stderrors::Nothing) = []
function sieve_results(estimate, filename, key, pair_to_var_id, sieve_stderrors)
    sieve_stderror = haskey(pair_to_var_id, (filename, key)) ? 
                            sieve_stderrors[pair_to_var_id[(filename, key)]] :
                            nothing
    return sieve_results_(estimate, sieve_stderror)
end

"""
If the influence curve is not saved, the report is already a summary (a NamedTuple)
"""
TMLE.summarize(s::NamedTuple) = s

function update_results!(results, tmlereport, phenotype_type, filename, key, pair_to_var_id, sieve_stderrors)
    s = summarize(tmlereport)
    phenotype = string(tmlereport.target_name)
    lwb, upb = s.confint

    summary = vcat(
        [phenotype, phenotype_type, s.query.name, filename, key,
        s.initial_estimate, s.estimate,
        s.stderror, s.pvalue, lwb, upb],
        sieve_results(s.estimate, filename, key, pair_to_var_id, sieve_stderrors)
    )

    push!(results, summary)
end

get_sieve_info(dirname_, sieve_filename_::Nothing) = (nothing, nothing)

function get_sieve_info(dirname_, sieve_filename_)
    sieve_file = jldopen(joinpath(dirname_, sieve_filename_))
    pair_to_var_id = Dict((f, r) => i for (i, (f, r)) in enumerate(sieve_file["SOURCEFILE_REPORTID_PAIRS"]))
    sieve_stderrors = sieve_file["STDERRORS"]
    close(sieve_file)
    return pair_to_var_id, sieve_stderrors
end

function build_summary(parsed_args)
    dirname_, prefix_ = splitdir(parsed_args["prefix"])
    dirname__ = dirname_ == "" ? "." : dirname_
    allhdf5files = filter(
            x -> startswith(x, prefix_) && endswith(x, ".hdf5"), 
            readdir(dirname__)
    )

    if parsed_args["sieve"]
        estimate_filenames = filter(x -> !occursin("sieve_variance", x) , allhdf5files)
        sieve_filename_ = filter(x -> occursin("sieve_variance", x) , allhdf5files)[1]
    else
        estimate_filenames = allhdf5files
        sieve_filename_ = nothing
    end
                            
    pair_to_var_id, sieve_stderrors = get_sieve_info(dirname_, sieve_filename_)

    results = init_results(sieve_filename_)
    for filename in estimate_filenames
        phenotype_type = occursin("Bool", filename) ? "BINARY" : "CONTINUOUS"
        jldopen(joinpath(dirname_, filename)) do io
            tmlereports = io["TMLEREPORTS"]
            for key in keys(tmlereports)
                tmlereport = tmlereports[key]
                update_results!(results, tmlereport, phenotype_type, filename, key, pair_to_var_id, sieve_stderrors)
            end
        end
    end

    CSV.write(parsed_args["out"], results)
end
