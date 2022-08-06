
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

function write_param_files(parsed_args,
                           param_files, 
                           binary_phenotypes,
                           continuous_phenotypes)
    batch_size = parsed_args["phenotype-batch-size"]
    binary_phenotypes_param_files = batched_param_files(param_files, binary_phenotypes, batch_size)
    continuous_phenotypes_param_files = batched_param_files(param_files, continuous_phenotypes, batch_size)

    for (i, param_file) in enumerate(binary_phenotypes_param_files)
        YAML.write_file(outpath(parsed_args, ".binary.parameter_$i.yaml"), param_file)
    end

    for (i, param_file) in enumerate(continuous_phenotypes_param_files)
        YAML.write_file(outpath(parsed_args, ".continuous.parameter_$i.yaml"), param_file)
    end 
end

"""
    validated_param_files(param_files, treatmentcols)

For now validation is only performed on treatment variables 
but could be xtended to other variables.
"""
function validated_param_files(param_files, treatmentcols)
    validated_param_files_ = Dict[]
    for param_file in param_files
        required_treatments = []
        update_treatments_list!(required_treatments, param_file["Treatments"])
        not_found = setdiff(required_treatments, treatmentcols)
        if length(not_found) == 0
            push!(validated_param_files_, param_file)
        else
            @warn string("Some treatment variables could not be read from the ",
                         "data files and associated parameter files will not ",
                         "be processed: ", join(not_found, ", "))
        end
    end
    return validated_param_files_
end

"""
    finalize_tmle_inputs(parsed_args)

Datasets are joined and rewritten to disk in the same order.
"""
function merge_and_write(parsed_args, genotypes, param_files)
    # Merge data and retrieve column names
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
    # Write data files
    for finalrole in ("binary-phenotypes", "continuous-phenotypes", "treatments", "confounders", "covariates")
        optionally_write_csv(outpath(parsed_args, ".$finalrole.csv"), final_dataset, columns[finalrole])
    end

    # Validate param files based on observed treatments
    param_files = validated_param_files(param_files, columns["treatments"])

    # Write param files
    write_param_files(
        parsed_args,
        param_files, 
        filter(!=("SAMPLE_ID"), columns["binary-phenotypes"]),
        filter(!=("SAMPLE_ID"), columns["continuous-phenotypes"])
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

update_treatments_list!(snp_list, other) = push!(snp_list, other)
function update_treatments_list!(snp_list, other::AbstractArray)
    for elem in other
        update_treatments_list!(snp_list, elem)
    end
end

isSNP(x) = occursin(r"^rs.*[0-9]+$", lowercase(x))

function snps_from_param_files(param_files)
    snps = String[]
    for param_file in param_files
        update_treatments_list!(snps, param_file["Treatments"])
    end
    return unique(filter(isSNP, snps))
end

function call_genotypes(probabilities::AbstractArray, variant_genotypes::AbstractVector, threshold::Real)
    n = size(probabilities, 2)
    t = Vector{Union{Int, Missing}}(missing, n)
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

all_snps_called(found_snps::Set, snp_list::AbstractVector) = Set(snp_list) == found_snps

"""
    genotypes_encoding(variant)

Since we are only considering bi-allelic variants, genotypes are encoded 
as the number of minor alleles.
"""
function genotypes_encoding(variant)
    minor = minor_allele(variant)
    all₁, _ = alleles(variant)
    if all₁ == minor
        return [2, 1, 0]
    else
        return [0, 1, 2]
    end
end

"""
    bgen_files(snps, bgen_prefix)

This function assumes the UK-Biobank structure
"""
function call_genotypes(bgen_prefix::String, snp_list::AbstractVector, threshold::Real)
    chr_dir_, prefix_ = splitdir(bgen_prefix)
    chr_dir = chr_dir_ == "" ? "." : chr_dir_
    genotypes = nothing
    found_snps = Set()
    for filename in readdir(chr_dir)
        all_snps_called(found_snps, snp_list) ? break : nothing
        if is_numbered_chromosome_file(filename, prefix_)
            bgenfile = read_bgen(joinpath(chr_dir_, filename))
            chr_genotypes = DataFrame(SAMPLE_ID=bgenfile.samples)
            for variant in BGEN.iterator(bgenfile)
                rsid_ = rsid(variant)
                if rsid_ ∈ snp_list
                    push!(found_snps, rsid_)
                    if n_alleles(variant) != 2
                        @warn("Skipping $rsid_, not bi-allelic")
                        continue
                    end
                    minor_allele_dosage!(bgenfile, variant)
                    variant_genotypes = genotypes_encoding(variant)
                    probabilities = probabilities!(bgenfile, variant)
                    chr_genotypes[!, rsid_] = call_genotypes(probabilities, variant_genotypes, threshold)
                end
            end
            genotypes = genotypes isa Nothing ? chr_genotypes :
                    innerjoin(genotypes, chr_genotypes, on=:SAMPLE_ID)
        end
    end
    return genotypes
end

function tmle_inputs(parsed_args)
    if parsed_args["%COMMAND%"] == "with-asb-trans"
        return 
    elseif parsed_args["%COMMAND%"] == "with-param-files"
        param_files = TargeneCore.load_param_files(parsed_args["with-param-files"]["param-prefix"])
        snp_list = TargeneCore.snps_from_param_files(param_files)
        genotypes = TargeneCore.call_genotypes(parsed_args["bgen-prefix"], snp_list, parsed_args["call-threshold"])
        merge_and_write(parsed_args, genotypes, param_files)
    else
        throw(ArgumentError("Unrecognized command."))
    end
end
