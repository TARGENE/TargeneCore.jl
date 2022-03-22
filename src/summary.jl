init_results() = DataFrame(
    PHENOTYPE = String[],
    PHENOTYPE_TYPE = String[],
    QUERYNAME = String[],
    INITIAL_ESTIMATE = Float64[],
    ESTIMATE = Float64[],
    STDERR = Float64[],
    PVAL = Float64[],
    LWB = Float64[],
    UPB = Float64[],
    SIEVE_STDERR = Union{Float64, Missing}[],
    SIEVE_PVAL = Union{Float64, Missing}[],
    SIEVE_LWB = Union{Float64, Missing}[],
    SIEVE_UPB = Union{Float64, Missing}[],
    )

sieve_results(estimate, sieve_stderror::Nothing) =
    missing, missing, missing, missing

function sieve_results(estimate, sieve_stderror)
    sieve_test = OneSampleZTest(estimate, sieve_stderror, 1)
    sieve_pval = pvalue(sieve_test)
    sieve_lwb, sieve_upb = confint(sieve_test)
    return sieve_stderror, sieve_pval, sieve_lwb, sieve_upb
end

function update_results!(results, qr, sieve_stderror, phenotype_type)
    bf = briefreport(qr)
    phenotype = string(qr.target_name)
    lwb, upb = bf.confint

    sieve_stderror, sieve_pval, sieve_lwb, sieve_upb = sieve_results(bf.estimate, sieve_stderror)

    push!(
        results, 
        [phenotype, phenotype_type, bf.query.name, 
        bf.initial_estimate, bf.estimate,
        bf.stderror, bf.pvalue, lwb, upb,
        sieve_stderror, sieve_pval, sieve_lwb, sieve_upb]
    )
end


function build_summary(parsed_args)
    dirname_, prefix_ = splitdir(parsed_args["prefix"])
    dirname__ = dirname_ == "" ? "." : dirname_
    allhdf5files = filter(
            x -> startswith(x, prefix_) && endswith(x, ".hdf5"), 
            readdir(dirname__)
    )
    sieve_filename_ = filter(x -> occursin("sieve_variance", x) , allhdf5files)[1]
    estimate_filenames = filter(x -> !occursin("sieve_variance", x) , allhdf5files)

    sieve_file = jldopen(joinpath(dirname_, sieve_filename_))
    pair_to_var_id = Dict((f, r) => i for (i, (f, r)) in enumerate(sieve_file["SOURCEFILE_REPORTID_PAIRS"]))
    sieve_stderrors = sieve_file["STDERRORS"]
    close(sieve_file)

    results = init_results()

    for filename in estimate_filenames
        phenotype_type = occursin("Bool", filename) ? "BINARY" : "CONTINUOUS"
        jldopen(joinpath(dirname_, filename)) do io
            qrs = haskey(io, "MACHINE") ? queryreports(io["MACHINE"]) : io["QUERYREPORTS"]
            for (qr_idx, qr) in enumerate(qrs)
                sieve_stderror = haskey(pair_to_var_id, (filename, qr_idx)) ? 
                    sieve_stderrors[pair_to_var_id[(filename, qr_idx)]] :
                    nothing
                update_results!(results, qr, sieve_stderror, phenotype_type)

            end
        end
    end

    CSV.write(parsed_args["out"], results)
end