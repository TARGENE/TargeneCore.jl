set_from_txt_file(filepath::AbstractString) = Set(open(readlines, filepath))

function make_dataset(genotypes_prefix,
    traits_file,
    pcs_file,
    variants_file;
    out="dataset.arrow",
    call_threshold=0.9,
    verbosity=0)

    # Read Traits
    verbosity > 0 && @info "Reading Traits."
    traits = read_csv_file(traits_file)

    # Read Confounders
    verbosity > 0 && @info "Reading Confounders."
    pcs = load_flash_pca_results(pcs_file)

    # Extract genotypes matching variants
    verbosity > 0 && @info "Reading Genotypes."
    genotypes = call_genotypes(
        genotypes_prefix, 
        set_from_txt_file(variants_file), 
        call_threshold
    )

    # Merge and write
    verbosity > 0 && @info "Merging and writing dataset."
    dataset = merge(traits, pcs, genotypes)
    Arrow.write(out, dataset)

    verbosity > 0 && @info "Done."
end
