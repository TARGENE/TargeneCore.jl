const CHR_REG = r"chr[1-9]+"

columnnames_no_sid(df) = filter(!=("SAMPLE_ID"), names(df))

const ROLE_MAPPING = [
    "genetic-confounders" => "confounders",
    "extra-confounders" => "confounders",
    "covariates" => "covariates",
    "extra-treatments" => "treatments"
]

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
            new_param_file["Targets"] = batch
            push!(new_param_files, new_param_file)
        end
    end
    return new_param_files
end

function write_param_files(outprefix,
                           param_files, 
                           binary_phenotypes,
                           continuous_phenotypes,
                           batch_size)
    binary_phenotypes_param_files = batched_param_files(param_files, binary_phenotypes, batch_size)
    continuous_phenotypes_param_files = batched_param_files(param_files, continuous_phenotypes, batch_size)

    for (i, param_file) in enumerate(binary_phenotypes_param_files)
        YAML.write_file(string(outprefix, ".binary.parameter_$i.yaml"), param_file)
    end

    for (i, param_file) in enumerate(continuous_phenotypes_param_files)
        YAML.write_file(string(outprefix, ".continuous.parameter_$i.yaml"), param_file)
    end 
end

"""
    validated_param_files(param_files, treatmentcols)

For now validation is only performed on treatment variables 
but could be xtended to other variables.
"""
function validated_param_files(param_files, treatments)
    validated_param_files_ = Dict[]
    for param_file in param_files
        required_treatments = []
        update_treatments_list!(required_treatments, param_file["Treatments"])
        not_found = setdiff(required_treatments, names(treatments))
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

function validated_param_files(eQTLs, bQTLs, treatments, param_prefix; positivity_constraint=0.01)
    # Build parameter files
    param_files = Dict[]
    for generator in treatment_tuple_generators(eQTLs, bQTLs, param_prefix)
        TargeneCore.build_asb_trans_param_files!(param_files, generator, treatments; positivity_constraint=positivity_constraint)
    end
    return param_files
end

"""
    finalize_tmle_inputs(parsed_args)

Datasets are joined and rewritten to disk in the same order.
"""
function build_final_dataset(parsed_args, genotypes)
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

    return final_dataset, columns

end

function write_tmle_inputs(outprefix, final_dataset, columns, param_files, batch_size)
    # Write data files
    for finalrole in ("binary-phenotypes", "continuous-phenotypes", "treatments", "confounders", "covariates")
        optionally_write_csv(string(outprefix, ".$finalrole.csv"), final_dataset, columns[finalrole])
    end
    write_param_files(
        outprefix,
        param_files, 
        filter(!=("SAMPLE_ID"), columns["binary-phenotypes"]),
        filter(!=("SAMPLE_ID"), columns["continuous-phenotypes"]),
        batch_size
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

update_treatments_list!(snp_list, other::Nothing) = nothing
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

function trans_actors(path)
    eQTLs = CSV.read(path, DataFrame, select=[:ID]).ID
    return unique(eQTLs)
end

function asb_snps(asb_prefix)
    bQTLs = []
    asbdir_, prefix_ = splitdir(asb_prefix)
    asbdir = asbdir_ == "" ? "." : asbdir_
    for file in readdir(asbdir)
        if occursin(prefix_, file)
            new_bQTLs = CSV.read(joinpath(asbdir_, file), DataFrame, select=[:ID, :CHROM, :isASB])
            for row in eachrow(new_bQTLs)
                if row.isASB === true && row.CHROM ∉ ("chrX", "chrY")
                    push!(bQTLs, row.ID)
                end
            end
        end
    end
    return unique(bQTLs)
end

function additive_interaction_settings(treatment_tuple, treatments::DataFrame)
    control_cases = []
    for t in treatment_tuple
        unique_vals = sort(unique(skipmissing(treatments[!, t])))
        push!(control_cases, collect(combinations(unique_vals, 2)))
    end
    return Iterators.product(control_cases...)
end

function satisfies_positivity(interaction_setting, freqs; positivity_constraint=0.01)
    for base_setting in Iterators.product(interaction_setting...)
        if !haskey(freqs, base_setting) || freqs[base_setting] < positivity_constraint
            return false
        end
    end
    return true
end

function frequency_table(treatments, treatment_tuple::AbstractVector)
    freqs = Dict()
    N = nrow(treatments)
    for (key, group) in pairs(groupby(treatments, treatment_tuple; skipmissing=true))
        freqs[values(key)] = nrow(group) / N
    end
    return freqs
end

function build_asb_trans_param_files!(param_files, treatment_tuple_generator, treatments; positivity_constraint=0.01)
    for treatment_tuple in treatment_tuple_generator
        treatment_tuple = collect(treatment_tuple)
        freqs = frequency_table(treatments, treatment_tuple)
        interaction_settings = additive_interaction_settings(treatment_tuple, treatments)

        param_file = Dict("Treatments" => treatment_tuple, "Parameters" => Dict[])
        for interaction_setting in interaction_settings
            if satisfies_positivity(interaction_setting, freqs; positivity_constraint=positivity_constraint)
                param = Dict{String, Any}("name" => "I")
                for var_index in eachindex(treatment_tuple)
                    control, case = interaction_setting[var_index]
                    treatment = treatment_tuple[var_index]
                    param["name"] = string(param["name"], "__", treatment, "_", control, "->", case)
                    param[treatment] = Dict("control" => control, "case" => case)
                end
                push!(param_file["Parameters"], param)
            end
        end

        # Push valid parameters
        if param_file["Parameters"] != Dict[]
            push!(param_files, param_file)
        end
    end
    return param_files
end


treatment_tuple_generators(eQTLs, bQTLs, param_prefix::Nothing) = (Iterators.product(eQTLs, bQTLs),)

function treatment_tuple_generators(eQTLs, bQTLs, param_prefix::String)
    param_files = TargeneCore.load_param_files(param_prefix)
    generators = []
    for param_file in param_files
        extra_treatments = []
        update_treatments_list!(extra_treatments, param_file["Treatments"])
        push!(generators, Iterators.product(eQTLs, bQTLs, [[x] for x in extra_treatments]...))
    end
    return generators
end

function tmle_inputs(parsed_args)
    batch_size = parsed_args["phenotype-batch-size"]
    outprefix = parsed_args["out-prefix"]
    call_threshold = parsed_args["call-threshold"]
    bgen_prefix = parsed_args["bgen-prefix"]
    positivity_constraint = parsed_args["positivity-constraint"]
    if parsed_args["%COMMAND%"] == "with-asb-trans"
        param_prefix = parsed_args["with-asb-trans"]["param-prefix"]
        # Retrive SNPs
        eQTLs = TargeneCore.trans_actors(parsed_args["with-asb-trans"]["trans-actors"])
        bQTLs = TargeneCore.asb_snps(parsed_args["with-asb-trans"]["asb-prefix"])
        snp_list = vcat(eQTLs, bQTLs)
        # Genotypes and final dataset
        genotypes = TargeneCore.call_genotypes(bgen_prefix, snp_list, call_threshold)
        final_dataset, columns = build_final_dataset(parsed_args, genotypes)
        # Parameter files
        param_files = validated_param_files(eQTLs, bQTLs, final_dataset[!, columns["treatments"]], param_prefix; 
                                            positivity_constraint=positivity_constraint)
        # write genotypes and queries
        write_tmle_inputs(outprefix, final_dataset, columns, param_files, batch_size)
    elseif parsed_args["%COMMAND%"] == "with-param-files"
        param_files = TargeneCore.load_param_files(parsed_args["with-param-files"]["param-prefix"])
        snp_list = TargeneCore.snps_from_param_files(param_files)
        # Genotypes and final dataset
        genotypes = TargeneCore.call_genotypes(bgen_prefix, snp_list, call_threshold)
        final_dataset, columns = build_final_dataset(parsed_args, genotypes)
        # Parameter files
        param_files = validated_param_files(param_files, final_dataset[!, columns["treatments"]])
        write_tmle_inputs(outprefix, final_dataset, columns, param_files, batch_size)
    else
        throw(ArgumentError("Unrecognized command."))
    end
end
