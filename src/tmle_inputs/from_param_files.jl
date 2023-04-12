
update_treatments_list!(snp_list, other::Nothing) = nothing

"""
    param_dicts_from_files(param_prefix)

Loads potentially not fully instantiated parameters files to 
dictionary representations.
"""
function param_dicts_from_files(param_prefix)
    param_files = Dict[]
    paramdir_, prefix = splitdir(param_prefix)
    paramdir = paramdir_ == "" ? "." : paramdir_
    for filename in readdir(paramdir)
        if startswith(filename, prefix)
            push!(param_files, YAML.load_file(joinpath(paramdir_, filename)))
        end
    end
    return param_files
end

function snps_and_variables_from_param_files(param_files, pcs, traits)
    snps = []
    variables = Dict(
        "PCs" => TargeneCore.pcnames(pcs),
        "extraT" => String[],
        "extraW" => String[],
        "extraC" => String[]
    )
    all_trait_names = Set(names(traits))
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
            if T ∈ all_trait_names
                unique!(push!(variables["extraT"], T))
            else
                unique!(push!(snps, T))
            end
        end
    end
    # Determine all potential `Y`s
    non_targets = Set(vcat("SAMPLE_ID", variables["extraT"], variables["extraW"], variables["extraC"]))
    variables["Y"] = targets_from_traits(traits, non_targets)

    return unique(snps), variables
end

MissingTraitVariablesError(variables) =
    ArgumentError(string("Requested trait variables: ", join(variables, ", "), " could not be found in the trait dataset."))

function maybe_throw_missing_traits(found, requested)
    if length(found) != length(requested)
        diff = setdiff(requested, found)
        throw(MissingTraitVariablesError(diff))
    end
end

AbsentAlleleError(variant, allele) =
    ArgumentError(string("Variant ", variant, "'s allele ", allele, " cannot be found in the data."))

"""
If the proposed allele is a String and cannot be found in the genotypes, 
if the reverse string allele can be found, then we replace it by this new 
representation. Otherwise, there is nothing that can be done.
"""
function try_fix_mismatch_and_report_success!(dict, key, actual_alleles)
    allele = dict[key]
    if !(allele ∈ actual_alleles)
        if allele isa Real
            return false
        elseif allele isa String
            length(allele) == 2 || throw("Only SNPs can be dealt with at the moment.")
            rev_allele = reverse(allele)
            if rev_allele ∈ actual_alleles
                dict[key] = rev_allele
            else
                return false
            end
        else
            throw(ArgumentError(string("Can't deal with allele ", allele, " of type: ", typeof(allele))))
        end
    end
    return true
end

"""
    adjust_treatment_encoding!(param_file, genotypes::DataFrame)

This function adjusts the `Parameters` section of a parameter file. 
    
Currently, there is one main issue to deal with. Since we are assuming 
the genotypes are unphased, if a variant's allele is represented 
as "AG" in the parameter file and as "GA" in the genotypes, there will be a mismatch 
between those representations that we want to resolve. This is done 
by modifying the parameter file's representation.

Furthermore if the parameter file's allele is simply not present in the 
genotypes, an error in thrown.
"""
function adjust_treatment_encoding!(param_file, genotypes::DataFrame)
    variants = Set(filter(x -> x != "SAMPLE_ID", string.(names(genotypes))))
    variant_alleles = Dict(v => unique(skipmissing(genotypes[!, v])) for v in variants)
    for param in param_file["Parameters"]
        for (key, val) in param
            if key ∈ variants
                adjust_treatment_encoding!(param, variant_alleles, key, val)
            end
        end
    end
end

"""
    adjust_treatment_encoding!(param, variant_alleles, key, val)

This adjusts the treatment encoding for parameters that don't have a case/control specification.
"""
function adjust_treatment_encoding!(param, variant_alleles, key, val)
    actual_alleles = variant_alleles[key]
    if !(val ∈ actual_alleles)
        s = try_fix_mismatch_and_report_success!(param, key, actual_alleles)
        s || throw(AbsentAlleleError(key, control_case_dict[c]))
    end
end

"""
    adjust_treatment_encoding!(param, variant_alleles, key, val)

This adjusts the treatment encoding for parameters that have a case/control specification.
"""
function adjust_treatment_encoding!(param, variant_alleles, key, control_case_dict::Dict)
    actual_alleles = variant_alleles[key]
    for c in ("case", "control")
        # If there is a mismatch between a case/control value in the 
        # proposed param file and the actual genotypes
        if !(control_case_dict[c] ∈ actual_alleles)
            s = try_fix_mismatch_and_report_success!(control_case_dict, c, actual_alleles)
            s || throw(AbsentAlleleError(key, control_case_dict[c]))
        end
    end
end


treatment_setting_from_param_dict(::Type{<:Union{ATE, IATE}}, param_dict::Dict, treatment_variables) = 
    [(param_dict[t]["control"], param_dict[t]["case"]) for t in treatment_variables]

treatment_setting_from_param_dict(::Type{CM}, param_dict::Dict, treatment_variables) = 
    [(param_dict[t],) for t in treatment_variables]

function param_dict_satisfies_positivity(param_dict, treatment_variables, freqs; positivity_constraint=0.)
    param_type = getfield(TMLE, Symbol(param_dict["name"]))
    treatment_setting = treatment_setting_from_param_dict(param_type, param_dict, treatment_variables)
    return satisfies_positivity(param_type, treatment_setting, freqs; positivity_constraint=positivity_constraint)
end

NoRemainingParamsError(positivity_constraint) = ArgumentError(string("No parameter passed the given positivity constraint: ", positivity_constraint))

function positivity_respecting_parameters(params_dict, data; positivity_constraint=0.)
    treatment_variables = params_dict["T"]
    freqs = TargeneCore.frequency_table(data, treatment_variables)
    filtered_param_dicts =  filter(
        param_dict -> param_dict_satisfies_positivity(
            param_dict, treatment_variables, freqs; 
            positivity_constraint=positivity_constraint
            ), 
        params_dict["Parameters"]
        )
    length(filtered_param_dicts) > 0 || throw(NoRemainingParamsError(positivity_constraint))
    return filtered_param_dicts
end

function adjusted_param_files!(param_dicts, variables, data, snp_list; positivity_constraint=0., batch_size=nothing)
    new_param_dicts = Dict[]
    for param_dict in param_dicts
        # If the genotypes encoding is a string representation make sure they match the actual genotypes
        TargeneCore.adjust_treatment_encoding!(param_dict, data[!, snp_list])
        # Update Parameters section
        param_dict["Parameters"] = positivity_respecting_parameters(param_dict, data; positivity_constraint=positivity_constraint)
        # Update `W` section with PCs
        if haskey(param_dict, "W")
            param_dict["W"] = unique(vcat(param_dict["W"], variables["PCs"]))
        else
            param_dict["W"] = variables["PCs"]
        end
        # If no target has been specified, use all available ones
        if !haskey(param_dict, "Y")
            add_batchified_param_files!(new_param_dicts, param_dict, variables["Y"], batch_size)
        else
            Y = intersect(param_dict["Y"], variables["Y"])
            maybe_throw_missing_traits(Y, param_dict["Y"])
            add_batchified_param_files!(new_param_dicts, param_dict, Y, batch_size)
        end
    end
    return new_param_dicts
end

MismatchedCaseControlEncodingError() = 
    ArgumentError(
                "All case/control values in the provided parameter 
                files should respect the same encoding. They should either 
                represent the number of minor alleles (i.e. 0, 1, 2) or use the genotypes string representation."
                )

function check_genotypes_encoding(val::Dict, type)
    if !(typeof(val["case"]) <: type && typeof(val["control"]) <: type)
        throw(MismatchedCaseControlEncodingError())
    end
end

check_genotypes_encoding(val::T, type) where T = 
    T <: type || throw(MismatchedCaseControlEncodingError())

get_genotype_encoding(val::Dict) = typeof(val["case"]) <: Real ? Real : String
get_genotype_encoding(val::T) where {T<:String} = String
get_genotype_encoding(val::T) where {T<:Real} = Real
get_genotype_encoding(val) = throw(ArgumentError(string("Genotype value(s):", val, " not allowed")))

"""
We enforce that all genotypes are encoded in the same way.
That is either with an integer representing the number of minor alleles 
or the string representation.
"""
function get_genotype_encoding(param_files, snps)
    genotypes_encoding = nothing
    for param_file in param_files
        for param in param_file["Parameters"]
            for (key, val) in param
                if key in snps
                    if genotypes_encoding === nothing 
                        genotypes_encoding = get_genotype_encoding(val)
                    else
                        check_genotypes_encoding(val, genotypes_encoding)
                    end
                end
            end
        end
    end
    return genotypes_encoding <: Real
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
    param_dicts = TargeneCore.param_dicts_from_files(param_prefix)

    # Genotypes and full data
    snp_list, variables = TargeneCore.snps_and_variables_from_param_files(param_dicts, pcs, traits)
    genotypes_asint = TargeneCore.get_genotype_encoding(param_dicts, snp_list)
    genotypes = TargeneCore.call_genotypes(bgen_prefix, snp_list, call_threshold; asint=genotypes_asint)
    data = TargeneCore.merge(traits, pcs, genotypes)

    # Parameter files
    param_dicts = TargeneCore.adjusted_param_files!(
        param_dicts, variables, data, snp_list; 
        positivity_constraint=positivity_constraint, batch_size=batch_size
    )

    # write data and parameter files
    write_tmle_inputs(outprefix, data, param_dicts)
end