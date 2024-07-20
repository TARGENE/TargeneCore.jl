function infer_method_from_config(config_file)
    if endswith(config_file, "jls")
        return tl_inputs_from_param_files
    elseif endswith(config_file, "yaml")
        config = YAML.load_file(config_file)
        if config["type"] == "Configuration" || config["type"] == "TMLE.Configuration"
            return tl_inputs_from_param_files
        else
            return allele_independent_estimands
        end
    else
        throw(ArgumentError("Unsupported file extension."))
    end
end

function estimation_inputs(config_file, genotypes_prefix, traits_file, pcs_file;
    outprefix="final",
    batchsize=nothing,
    call_threshold=0.9,
    positivity_constraint=0.01,
    verbosity=0
    )
    estimation_inputs_method = infer_method_from_config(config_file)
    estimation_inputs_method(config_file, genotypes_prefix, traits_file, pcs_file;
        outprefix=outprefix,
        batchsize=batchsize,
        call_threshold=call_threshold,
        positivity_constraint=positivity_constraint,
        verbosity=verbosity
    )
    return 0
end