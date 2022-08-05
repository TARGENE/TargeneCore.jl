
columnnames_no_sid(df) = filter(!=("SAMPLE_ID"), names(df))

const ROLE_MAPPING = [
    "genetic-confounders" => "confounders",
    "extra-confounders" => "confounders",
    "covariates" => "covariates",
    "extra-treatments" => "treatments"
]

outpath(parsed_args, suffix) = string(parsed_args["out-prefix"], suffix)

function optionally_write_csv(path, final_dataset, columns)
    if columns != ["SAMPLE_ID"]
        CSV.write(path, final_dataset[!, columns])
    end
end

function update_with_extra_treatment!(param_file, extra_treatments)
    append!(param_file["treatments"], extra_treatments)
    for param in param_file["Parameters"]

    end
end

function write_phenotype_batch(phenotypes, batch, batch_start, batch_end)
    batch_file = batch_filename(batch)
    open(batch_file, "w") do io
        for i in batch_start:batch_end
            println(io, phenotypes[i])
        end
    end
end

batch_filename(batch) = string("phenotypes_batch_", batch, ".csv")

"""

If no batch-size is provided a batch consists of all phenotypes
"""
batched_param_files(param_files, phenotypes, batch_size::Nothing) = param_files

function batched_param_files(param_files, phenotypes, batch_size::Int)
    new_param_files = Dict[]
    for param_file in param_files
        for batch in Iterators.partition(phenotypes, batch_size)
            new_param_file = copy(param_file)
            new_param_file["Phenotypes"] = batch
            push!(new_param_files, new_param_file)
        end
    end
    return new_param_files
end

function write_param_files(param_files, 
                           rsid_to_minor_major,
                           extra_treatments,
                           batch_size,
                           binary_phenotypes,
                           continuous_phenotypes)

    binary_phenotypes_param_files = batched_param_files(param_files, binary_phenotypes, batch_size)
    continuous_phenotypes_param_files = batched_param_files(param_files, continuous_phenotypes, batch_size)
    for param_file in vcat(binary_phenotypes_param_files, continuous_phenotypes_param_files)
        update_treatment_section(param_file, extra_treatments)
    end

end
"""
    finalize_tmle_inputs(parsed_args)

Datasets are joined and rewritten to disk in the same order.
"""
function merge_and_write(parsed_args, genotypes, rsid_to_minor_major, param_files)
    binary_phenotypes = CSV.read(parsed_args["binary-phenotypes"], DataFrame)
    continuous_phenotypes = CSV.read(parsed_args["continuous-phenotypes"], DataFrame)
    columns = Dict(
        "binary-phenotypes" => names(binary_phenotypes),
        "continuous-phenotypes" => names(continuous_phenotypes),
        "covariates" => ["SAMPLE_ID"],
        "treatments" => names(genotypes),
        "confounders" => ["SAMPLE_ID"]
        )
    final_dataset = innerjoin(
        innerjoin(binary_phenotypes, continuous_phenotypes, on=:SAMPLE_ID), 
        genotypes,
        on=:SAMPLE_ID
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
    optionally_write_csv(outpath(parsed_args, ".binary-phenotypes.csv"), final_dataset, columns["binary-phenotypes"])
    optionally_write_csv(outpath(parsed_args, ".continuous-phenotypes.csv"), final_dataset, columns["continuous-phenotypes"])
    # Write treatments
    optionally_write_csv(outpath(parsed_args, ".treatments.csv"), final_dataset, columns["treatments"])
    # Write confounders
    optionally_write_csv(outpath(parsed_args, ".confounders.csv"), final_dataset, columns["confounders"])
    # Write covariates
    optionally_write_csv(outpath(parsed_args, ".covariates.csv"), final_dataset, columns["covariates"])
    # Write param files
    write_param_files(
        param_files, 
        rsid_to_minor_major,
        filter(x -> x ∉ names(genotypes), columns["treatments"]),
        parsed_args["phenotype-batch-size"],
        columns["binary-phenotypes"],
        columns["continuous-phenotypes"]
        )
end

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

update_snp_list!(snp_list, other) = push!(snp_list, other)
function update_snp_list!(snp_list, other::AbstractArray)
    for elem in other
        update_snp_list!(snp_list, elem)
    end
end

function snps_from_param_files(param_files)
    snps = String[]
    for param_file in param_files
        update_snp_list!(snps, param_file["Treatments"])
    end
    return unique(snps)
end

function call_genotypes(probabilities::AbstractArray, variant_genotypes::AbstractVector, threshold::Real)
    n = size(probabilities, 2)
    t = Vector{Union{String, Missing}}(missing, n)
    for i in 1:n
        # If no allele has been annotated with sufficient confidence
        # the sample is declared as missing for this variant
        genotype_index = findfirst(x -> x >= threshold, probabilities[:, i])
        genotype_index isa Nothing || (t[i] = variant_genotypes[genotype_index])
    end
    return t
end

function is_numbered_chromosome_file(filename, prefix)
    if occursin(prefix, filename) && endswith(filename, "bgen")
        regexp_match = match(CHR_REG, filename)
        if regexp_match !== nothing
            return true
        end
    end
    return false
end

function read_bgen(filepath)
    sample_filepath = string(filepath[1:end-4], "sample")
    idx_filepath = string(filepath, ".bgi")
    return Bgen(filepath, sample_path=sample_filepath, idx_path=idx_filepath)
end

all_snps_called(genotypes::Nothing, sns_list::AbstractVector) = false
function all_snps_called(genotypes::DataFrame, snp_list::AbstractVector)
    if size(genotypes, 2) == size(snp_list, 1) + 1
        return true
    end
    return false
end

"""
    bgen_files(snps, bgen_prefix)

This assumes the UK-Biobank structure
"""
function call_genotypes(bgen_prefix::String, snp_list::AbstractVector, threshold::Real)
    chr_dir_, prefix_ = splitdir(bgen_prefix)
    chr_dir = chr_dir_ == "" ? "." : chr_dir_
    rsid_to_minor_major = Dict()
    genotypes = nothing
    for filename in readdir(chr_dir)
        all_snps_called(genotypes, snp_list) ? break : nothing
        if is_numbered_chromosome_file(filename, prefix_)
            bgenfile = read_bgen(joinpath(chr_dir_, filename))
            chr_genotypes = DataFrame(SAMPLE_ID=bgenfile.samples)
            for variant in BGEN.iterator(bgenfile)
                rsid_ = rsid(variant)
                if rsid_ ∈ snp_list
                    if n_alleles(variant) != 2
                        @warn("Skipping $rsid_, not bi-allelic")
                        continue
                    end
                    minor_allele_dosage!(bgenfile, variant)
                    major = major_allele(variant)
                    minor = minor_allele(variant)
                    all₁, all₂ = alleles(variant)
                    variant_genotypes = [all₁*all₁, major*minor, all₂*all₂]
                    rsid_to_minor_major[rsid_] = (minor=minor, major=major)
                    probabilities = probabilities!(bgenfile, variant)
                    chr_genotypes[!, rsid_] = call_genotypes(probabilities, variant_genotypes, threshold)
                end
            end
            genotypes = genotypes isa Nothing ? chr_genotypes :
                    innerjoin(genotypes, chr_genotypes, on=:SAMPLE_ID)
        end
    end
    return genotypes, rsid_to_minor_major
end

function tmle_inputs(parsed_args)
    if parsed_args["%COMMAND%"] == "with-asb-trans"
        return 
    elseif parsed_args["%COMMAND%"] == "with-param-files"
        param_files = TargeneCore.load_param_files(parsed_args["with-param-files"]["param-prefix"])
        snp_list = TargeneCore.snps_from_param_files(param_files)
        genotypes, rsid_to_minor_major = TargeneCore.call_genotypes(parsed_args["bgen-prefix"], snp_list, parsed_args["call-threshold"])
        merge_and_write(parsed_args, genotypes, rsid_to_minor_major, param_files)
    else
        throw(ArgumentError("Unrecognized command."))
    end
end