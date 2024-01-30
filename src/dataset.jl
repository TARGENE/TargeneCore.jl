function generate_dataset(parsed_args)
    verbosity = parsed_args["verbosity"]

    # Read Traits
    verbosity > 0 && @info "Reading Traits."
    traits = read_csv_file(parsed_args["traits"])

    # Read Confounders
    verbosity > 0 && @info "Reading Confounders."
    pcs = read_csv_file(parsed_args["confounders"])

    # Extract genotypes matching variants
    verbosity > 0 && @info "Reading Genotypes."
    genotypes = call_genotypes(
        parsed_args["bgen-prefix"], 
        set_from_txt_file(parsed_args["variants"]), 
        parsed_args["call-threshold"]
    )

    # Merge and write
    verbosity > 0 && @info "Merging and writing dataset."
    dataset = merge(traits, pcs, genotypes)
    Arrow.write(parsed_args["out"], dataset)

    verbosity > 0 && @info "Done."
end
