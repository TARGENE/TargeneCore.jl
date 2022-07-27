
columnnames_no_sid(df) = filter(!=("SAMPLE_ID"), names(df))

prepend_sample_id(vec) = vcat("SAMPLE_ID", vec)

const ROLE_MAPPING = [
    "genetic-confounders" => "confounders",
    "extra-confounders" => "confounders",
    "covariates" => "covariates",
    "genotypes" => "treatments",
    "extra-treatments" => "treatments"
]

"""
    finalize_tmle_inputs(parsed_args)

Datasets are joined and rewritten to disk in the same order.
"""
function finalize_tmle_inputs(parsed_args)
    final_dataset = CSV.read(parsed_args["phenotypes"], DataFrame)
    columns = Dict(
        "phenotypes" => columnnames_no_sid(final_dataset),
        "covariates" => String[],
        "treatments" => String[],
        "confounders" => String[]
        )
    for (key, role) in ROLE_MAPPING
        if parsed_args[key] !== nothing
            new_data = CSV.read(parsed_args[key], DataFrame)
            append!(columns[role], columnnames_no_sid(new_data))
            final_dataset = innerjoin(
                final_dataset,
                new_data,
                on=:SAMPLE_ID
            )
        end
    end

    # Write targets
    CSV.write(
        string(parsed_args["out-prefix"], ".phenotypes.csv"), 
        final_dataset[!, prepend_sample_id(columns["phenotypes"])]
    )
    # Write treatments
    CSV.write(
        string(parsed_args["out-prefix"], ".treatments.csv"), 
        final_dataset[!, prepend_sample_id(columns["treatments"])]
    )
    # Write confounders
    CSV.write(
        string(parsed_args["out-prefix"], ".confounders.csv"), 
        final_dataset[!, prepend_sample_id(columns["confounders"])]
    )
    # Write covariates
    if columns["covariates"] != String[]
        CSV.write(
            string(parsed_args["out-prefix"], ".covariates.csv"), 
            final_dataset[!, prepend_sample_id(columns["covariates"])]
        )
    end
end