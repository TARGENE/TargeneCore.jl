
update_treatments_list!(snp_list, other::Nothing) = nothing
function load_param_files(param_prefix)
    param_files = Dict[]
    paramdir_, prefix = splitdir(param_prefix)
    paramdir = paramdir_ == "" ? "." : paramdir_
    for filename in readdir(paramdir)
        if occursin(prefix, filename)
            push!(param_files, YAML.load_file(joinpath(paramdir_, filename)))
        end
    end
    return param_files
end

isSNP(x) = occursin(r"^rs.*[0-9]+$", lowercase(x))

function snps_and_variables_from_param_files(param_files, pcs, traits)
    snps = []
    variables = Dict(
        "PCs" => TargeneCore.pcnames(pcs),
        "extraT" => String[],
        "extraW" => String[],
        "extraC" => String[]
    )
    for param_file in param_files
        # Parse `C` section
        if haskey(param_file, "C")
            union!(variables["extraC"], string.(param_file["C"]))
        end
        # Parse `W` section
        if haskey(param_file, "W")
            union!(variables["extraW"], string.(param_file["W"]))
        end
        # Parse `T` section
        for T in string.(param_file["T"])
            if TargeneCore.isSNP(T)
                unique!(push!(snps, T))
            else
                unique!(push!(variables["extraT"], T))
            end
        end
    end
    # Determine all potential `Y`s
    non_targets = Set(vcat("SAMPLE_ID", variables["extraT"], variables["extraW"], variables["extraC"]))
    variables["Y"] = targets_from_traits(traits, non_targets)

    return unique(snps), variables
end

function adjusted_param_files(param_files, variables; positivity_constraint=0., batch_size=nothing)
    new_param_files = Dict[]
    for param_file in param_files
        (haskey(param_file, "Parameters") && param_file["Parameters"] != Dict[]) ? nothing :
            throw(ArgumentError("One of the parameter file is incomplete, at least one parameter should be provided in the 'Parameters' section."))
        # Update `W` section with PCs
        if haskey(param_file, "W")
            param_file["W"] = unique(vcat(param_file["W"], variables["PCs"]))
        else
            param_file["W"] = variables["PCs"]
        end

        # If no target jas been specified, use all available ones
        if !haskey(param_file, "Y")
            add_batchified_param_files!(new_param_files, param_file, variables["Y"], batch_size)
        else
            Y = intersect(param_file["Y"], variables["Y"]) 
            add_batchified_param_files!(new_param_files, param_file, Y, batch_size)
        end
    end
    return new_param_files
end

function tmle_inputs_from_param_files(parsed_args)
    # Read parsed_args
    batch_size = parsed_args["phenotype-batch-size"]
    outprefix = parsed_args["out-prefix"]
    call_threshold = parsed_args["call-threshold"]
    bgen_prefix = parsed_args["bgen-prefix"]
    positivity_constraint = parsed_args["positivity-constraint"]
    traits = TargeneCore.read_data(parsed_args["traits"])
    pcs = TargeneCore.read_data(parsed_args["pcs"])
    param_prefix = parsed_args["from-param-files"]["param-prefix"]

    # Load initial Parameter files
    param_files = TargeneCore.load_param_files(param_prefix)

    # Genotypes and full data
    snp_list, variables = TargeneCore.snps_and_variables_from_param_files(param_files, pcs, traits)
    genotypes = TargeneCore.call_genotypes(bgen_prefix, snp_list, call_threshold)
    data = TargeneCore.merge(traits, pcs, genotypes)

    # Parameter files
    param_files = TargeneCore.adjusted_param_files(param_files, variables; positivity_constraint=positivity_constraint, batch_size=batch_size)

    # write data and parameter files
    write_tmle_inputs(outprefix, data, param_files)
end