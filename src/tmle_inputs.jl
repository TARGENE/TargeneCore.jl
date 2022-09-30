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
    new_param_files = Pair{String, Dict}[]
    for (filename, param_file) in param_files
        for (batch_id, batch) in enumerate(Iterators.partition(phenotypes, batch_size))
            new_param_file = copy(param_file)
            new_param_file["Targets"] = batch
            new_filename = string(replace(filename, ".yaml" => ""), ".batch_", batch_id, ".yaml")
            push!(new_param_files, new_filename => new_param_file)
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

    for (filename, param_file) in binary_phenotypes_param_files
        YAML.write_file(string(outprefix, ".binary.", filename), param_file)
    end

    for (filename, param_file) in continuous_phenotypes_param_files
        YAML.write_file(string(outprefix, ".continuous.", filename), param_file)
    end 
end

"""
    validated_param_files(param_files, treatmentcols)

For now validation is only performed on treatment variables 
but could be xtended to other variables.
"""
function validated_param_files(param_files, treatments)
    validated_param_files_ = Pair{String, Dict}[]
    for (filename, param_file) in param_files
        required_treatments = []
        update_treatments_list!(required_treatments, param_file["Treatments"])
        not_found = setdiff(required_treatments, names(treatments))
        if length(not_found) == 0
            push!(validated_param_files_, filename => param_file)
        else
            @warn string("Some treatment variables could not be read from the ",
                         "data files and associated parameter files will not ",
                         "be processed: ", join(not_found, ", "))
        end
    end
    return validated_param_files_
end


"""
    read_data(filepath)

The SAMPLE_ID column should be read as a String.
"""
read_data(filepath) = CSV.read(filepath, DataFrame, types=Dict(:SAMPLE_ID => String))


function build_final_dataset(parsed_args, genotypes, envs)
    # Merge data and retrieve column names
    traits = read_data(parsed_args["traits"])
    pcs = read_data(parsed_args["pcs"])
    non_targets = []

    # Confounders
    variables = Dict(
        "W" => filter(!==("SAMPLE_ID"), names(pcs)),
        "Y_binary" => [],
        "Y_continuous" => []
    )
    if haskey(parsed_args, "extra-confounders")
        extra_confounders = CSV.read(parsed_args["extra-confounders"], DataFrame, header=false)[!, 1]
        append!(variables["W"], extra_confounders)
        append!(non_targets, variables["W"])
    end
    
    # Extra covariates
    if haskey(parsed_args, "extra-covariates")
        extra_covariates = CSV.read(parsed_args["extra-covariates"], DataFrame, header=false)[!, 1]
        columns["C"] = extra_covariates
        append!(non_targets, extra_covariates)
    end

    # Extra covariates
    if envs !== nothing
        append!(non_targets, envs.ID)
    end

    for colname in names(traits)
        if colname ∉ non_targets
            col = traits[!, colname]
            # Targets must be either continuous or binary
            if Set(unique(skipmissing(x))) == Set([0, 1])
                push!(variables["Y_binary"], colname)
            else
                @assert nonmissingtype(eltype(col)) <: Real 
                    "Variable $colname is neither binary nor continuous "*
                    "but is tried to be used as a target. Make sure it is "
                    "correctly declared as a treatment, extra-covariate or "*
                    "extra-confounder."
                push!(variables["Y_continuous"], colname)
            end
        end
    end

    data = innerjoin(
        innerjoin(traits, pcs, on="SAMPLE_ID"),
        genotypes,
        on="SAMPLE_ID"
    )

    return data, variables

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
    param_files = Pair{String, Dict}[]
    paramdir_, prefix = splitdir(param_prefix)
    paramdir = paramdir_ == "" ? "." : paramdir_
    for filename in readdir(paramdir)
        if occursin(prefix, filename)
            push!(param_files, filename => YAML.load_file(joinpath(paramdir_, filename)))
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
    for (filename, param_file) in param_files
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

function validated_param_files(param_files, data, variables; positivity_constraint=0.01, batch_size=1)
    new_param_files = []
    for param_file in param_files
        treatment_tuple = param_file["T"]
        freqs = frequency_table(data, treatment_tuple)
        interaction_settings = additive_interaction_settings(treatment_tuple, data)
        param_file["Parameters"] = Dict[]
        for interaction_setting in interaction_settings
            if satisfies_positivity(interaction_setting, freqs; positivity_constraint=positivity_constraint)
                param = Dict{String, Any}("name" => "IATE")
                for var_index in eachindex(treatment_tuple)
                    control, case = interaction_setting[var_index]
                    treatment = treatment_tuple[var_index]
                    param[treatment] = Dict("control" => control, "case" => case)
                end
                push!(param_file["Parameters"], param)
            end
        end

        # Push valid parameters
        if param_file["Parameters"] != Dict[]
            param_file["W"] = variables["W"]
            param_file["C"] = variables["C"]
            for targettype in ("Y_continuous", "Y_binary")
                if batch_size !== nothing
                    for batch in Iterators.partition(variables[targettype], batch_size)
                        batched_param_file = copy(param_file)
                        batched_param_file["Y"] = batch
                        push!(new_param_files, batched_param_file)
                    end
                else
                    batched_param_file = copy(param_file)
                    batched_param_file["Y"] = variables[targettype]
                    push!(new_param_files, batched_param_file)
                end
            end
        end
    end
    return new_param_files
end

read_snps_from_csv(path::Nothing) = nothing
read_snps_from_csv(path::String) = unique(CSV.read(path, DataFrame; select=[:ID, :CHR]), :ID)

read_envs_from_file(path::Nothing) = nothing
read_envs_from_file(path::String) = CSV.read(path, DataFrame; header=["ENV"])

trans_actors_from_prefix(trans_actors_prefix::Nothing) = nothing
function trans_actors_from_prefix(trans_actors_prefix::String)
    dir, prefix = splitdir(trans_actors_prefix)
    dir_ = dir == "" ? "." : dir
    trans_actors = DataFrame[]
    for filename in readdir(dir_)
        if startswith(filename, prefix)
            push!(trans_actors, read_snps_from_csv(joinpath(dir, filename)))
        end
    end
    return trans_actors
end

actors_error() = throw(ArgumentError("At least two of (bqtls, env, trans-actors) should be specified"))

treatments_from_actors(::Nothing, ::Nothing, _) = actors_error()
treatments_from_actors(_, ::Nothing, ::Nothing) = actors_error()
treatments_from_actors(::Nothing, _, ::Nothing) = actors_error()

function treatments_from_actors(bqtl_file, env_file, trans_actors_prefix)
    bqtls = read_snps_from_csv(bqtl_file)
    envs  = read_envs_from_file(env_file)
    transactors = trans_actors_from_prefix(trans_actors_prefix)
    bqtls, transactors, envs
end

all_snps(bqtls::DataFrame, transactors::Nothing) = bqtls.ID
all_snps(bqtls::Nothing, transactors::Vector{DataFrame}) = vcat((x.ID for x in transactors)...)
all_snps(bqtls::DataFrame, transactors::Vector{DataFrame}) = vcat(bqtls.ID, (x.ID for x in transactors)...)


combine_trans_actors(trans_actors::Vector{DataFrame}, envs::DataFrame, order) = combinations([trans_actors..., envs], order)
combine_trans_actors(trans_actors::Vector{DataFrame}, envs::Nothing, order) = combinations(trans_actors, order)
combine_trans_actors(trans_actors::Nothing, envs::DataFrame, order) = [[envs]]

function combine_by_bqtl(bqtls::DataFrame, trans_actors::Union{Vector{DataFrame}, Nothing}, envs::Union{DataFrame, Nothing}, order::Int)
    causal_model_configurations = Dict[]
    for comb in combine_trans_actors(trans_actors, envs, order - 1)
        comb_interactions = crossjoin(bqtls, comb..., makeunique=true)
        for row in eachrow(comb_interactions)
            push!(
                causal_model_configurations,
                Dict(
                    "T" => vcat(row["ID"], [row[string("ID_", i)] for i in 1:order-1])
                )
            )
        end
    end
    return causal_model_configurations
end


function tmle_inputs(parsed_args)
    batch_size = parsed_args["phenotype-batch-size"]
    outprefix = parsed_args["out-prefix"]
    call_threshold = parsed_args["call-threshold"]
    bgen_prefix = parsed_args["bgen-prefix"]
    positivity_constraint = parsed_args["positivity-constraint"]
    if parsed_args["%COMMAND%"] == "with-actors"
        param_prefix = parsed_args["with-actors"]["param-prefix"]
        # Retrieve SNPs and environmental treatments
        bqtls, transactors, envs = treatments_from_actors(
            parsed_args["with-actors"]["bqtl"], 
            parsed_args["with-actors"]["env"], 
            parsed_args["with-actors"]["trans-actors"]
        )
        param_files = Dict[]
        for order in parsed_args["with-actors"]["orders"]
            append!(
                param_files,
                combine_by_bqtl(bqtls, trans_actors, envs, order)
            )
        end

        snp_list = all_snps(bqtls, transactors)
        # Genotypes and final dataset
        genotypes = TargeneCore.call_genotypes(bgen_prefix, snp_list, call_threshold)

        final_dataset, variables = build_final_dataset(parsed_args, genotypes, envs)

        # Parameter files
        param_files = validated_param_files(param_files, final_dataset; positivity_constraint=positivity_constraint)
        # write genotypes and queries
        write_tmle_inputs(outprefix, final_dataset, columns, param_files, batch_size)
    elseif parsed_args["%COMMAND%"] == "with-param-files"
        param_files = TargeneCore.load_param_files(parsed_args["with-param-files"]["param-prefix"])
        snp_list = TargeneCore.snps_from_param_files(param_files)
        # Genotypes and final dataset
        genotypes = TargeneCore.call_genotypes(bgen_prefix, snp_list, call_threshold)
        final_dataset, columns = TargeneCore.build_final_dataset(parsed_args, genotypes)
        # Parameter files
        param_files = TargeneCore.validated_param_files(param_files, final_dataset[!, columns["treatments"]])
        write_tmle_inputs(outprefix, final_dataset, columns, param_files, batch_size)
    else
        throw(ArgumentError("Unrecognized command."))
    end
end
