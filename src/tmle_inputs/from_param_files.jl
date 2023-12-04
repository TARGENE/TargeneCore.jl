MissingVariableError(variable) =
    ArgumentError(string("Requested variable: ", variable, " cannot be found in the trait dataset."))

NotSNPAndWontFixError(variant, allele) = 
    ArgumentError(string(variant, "'s allele: ", allele, " cannot be found the genotyping dataset."))

AbsentAlleleError(variant, allele) =
    ArgumentError(string("Variant ", variant, "'s allele ", allele, " cannot be found in the genotyping dataset."))

MismatchedCaseControlEncodingError() = 
    ArgumentError(
                "All case/control values in the provided parameter 
                files should respect the same encoding. They should either 
                represent the number of minor alleles (i.e. 0, 1, 2) or use the genotypes string representation."
                )

NoRemainingParamsError(positivity_constraint) = ArgumentError(string("No parameter passed the given positivity constraint: ", positivity_constraint))


function check_genotypes_encoding(val::NamedTuple, type)
    if !(typeof(val.case) <: type && typeof(val.control) <: type)
        throw(MismatchedCaseControlEncodingError())
    end
end

check_genotypes_encoding(val::T, type) where T = 
    T <: type || throw(MismatchedCaseControlEncodingError())


function get_variables(estimands, traits, pcs)
    genetic_variants = Set{Symbol}()
    others = Set{Symbol}()
    pcs = Set{Symbol}(filter(x -> x != :SAMPLE_ID, propertynames(pcs)))
    alltraits = Set{Symbol}(filter(x -> x != :SAMPLE_ID, propertynames(traits)))
    for Ψ in estimands
        treatments = keys(Ψ.treatment_values)
        confounders = Iterators.flatten((Tconf for Tconf ∈ Ψ.treatment_confounders))
        push!(
            others, 
            Ψ.outcome_extra_covariates..., 
            confounders..., 
            treatments...
        )
        for T in treatments
            T ∉ alltraits && push!(genetic_variants, T)
        end
    end

    outcomes = setdiff(alltraits, others)

    return (genetic_variants=genetic_variants, outcomes=outcomes, pcs=pcs)
end

fix_mismatch(variant, allele::T, actual_alleles) where T = 
    throw(ArgumentError(string("Can't deal with ", variant, "'s allele ", allele, " of type: ", T)))

function fix_mismatch(variant, allele::String, actual_alleles)
    length(allele) == 2 || throw(NotSNPAndWontFixError(variant, allele))
    rev_allele = reverse(allele)
    rev_allele ∈ actual_alleles && return rev_allele
    throw(AbsentAlleleError(variant, allele))
end

function adjust_treatment_encoding(treatment_setting, variants_alleles, variant)
    variant_alleles = variants_alleles[variant]
    if treatment_setting ∉ variant_alleles 
        return fix_mismatch(variant, treatment_setting, variant_alleles)
    end
    return treatment_setting
end


function adjust_treatment_encoding(treatment_setting::NamedTuple{names, }, variants_alleles, variant) where names
    variant_alleles = variants_alleles[variant]
    case_control_values = []
    for c in names
        # If there is a mismatch between a case/control value in the 
        # proposed param file and the actual genotypes
        if treatment_setting[c] ∉ variant_alleles
            new_setting = fix_mismatch(variant, treatment_setting[c], variant_alleles)
            push!(case_control_values, new_setting)
        else
            push!(case_control_values, treatment_setting[c])
        end
    end
    return NamedTuple{names}(case_control_values)
end

"""
    adjust_parameter_sections(Ψ::T, variants_alleles, pcs) where T<:TMLE.Estimand

This function adjusts the treatment field of a `TMLE.Estimand`. 
    
Currently, there is one main issue to deal with. Since we are assuming 
the genotypes are unphased, if a variant's allele is represented 
as "AG" in the parameter file and as "GA" in the genotypes, there will be a mismatch 
between those representations that we want to resolve.

Furthermore if the parameter file's allele is simply not present in the 
genotypes, an error in thrown.
"""
function adjust_parameter_sections(Ψ::T, variants_alleles, pcs) where T<:TMLE.Estimand
    treatment_variables = keys(Ψ.treatment_values)
    treatment_values = []
    treatment_confounders = []
    for treatment in treatment_variables
        treatment_setting = Ψ.treatment_values[treatment]
        # Update treatment's values encoding
        if treatment ∈ keys(variants_alleles)
            push!(
                treatment_values,
                adjust_treatment_encoding(treatment_setting, variants_alleles, treatment)
            )
        else
            push!(treatment_values, treatment_setting)
        end
        # Update treatment's confoudner's values with PCs
        push!(treatment_confounders, vcat(Ψ.treatment_confounders[treatment]..., pcs...))
    end
    treatments = NamedTuple{treatment_variables}(treatment_values)
    confounders = NamedTuple{treatment_variables}(treatment_confounders)

    return T(outcome=Ψ.outcome, treatment_values=treatments, treatment_confounders=confounders, outcome_extra_covariates=Ψ.outcome_extra_covariates)
end

    
function append_from_valid_parameters!(
    parameters::Vector{<:TMLE.Estimand},
    frequency_tables::Dict,
    Ψ::T, 
    data,
    variants_alleles,
    variables;
    positivity_constraint=0.
    ) where T
    # Update treatment's and confounders's sections of Ψ
    Ψ = adjust_parameter_sections(Ψ, variants_alleles, variables.pcs)
    # Update frequency tables with current treatments
    treatments = sorted_treatment_names(Ψ)
    if !haskey(frequency_tables, treatments)
        frequency_tables[treatments] = TargeneCore.frequency_table(data, collect(treatments))
    end
    # Check if parameter satisfies positivity
    satisfies_positivity(Ψ, frequency_tables[treatments]; 
        positivity_constraint=positivity_constraint) || return
    # Expand wildcard to all outcomes
    if Ψ.outcome === :ALL
        update_parameters_from_outcomes!(parameters, Ψ, variables.outcomes)
    else
        # Ψ.target || MissingVariableError(variable)
        push!(parameters, Ψ)
    end
end

function adjusted_parameters(estimands, variables, data; positivity_constraint=0.)
    final_estimands = TMLE.Estimand[]
    variants_alleles = Dict(v => Set(unique(skipmissing(data[!, v]))) for v in variables.genetic_variants)
    freqency_tables = Dict()
    for Ψ in estimands
        # If the genotypes encoding is a string representation make sure they match the actual genotypes
        append_from_valid_parameters!(final_estimands, freqency_tables, Ψ, data, variants_alleles, variables; positivity_constraint=positivity_constraint)
    end

    length(final_estimands) > 0 || throw(NoRemainingParamsError(positivity_constraint))
    ordered_estimands = groups_ordering(final_estimands; brute_force=false, do_shuffle=false)
    
    return ordered_estimands
end


function tmle_inputs_from_param_files(parsed_args)
    # Read parsed_args
    batch_size = parsed_args["batch-size"]
    outprefix = parsed_args["out-prefix"]
    call_threshold = parsed_args["call-threshold"]
    bgen_prefix = parsed_args["bgen-prefix"]
    positivity_constraint = parsed_args["positivity-constraint"]
    traits = TargeneCore.read_data(parsed_args["traits"])
    pcs = TargeneCore.read_data(parsed_args["pcs"])
    estimands_config_file = parsed_args["from-param-file"]["paramfile"]

    # Load initial Parameter files
    estimands = TMLE.read_yaml(estimands_config_file).estimands

    # Genotypes and full data
    variables = TargeneCore.get_variables(estimands, traits, pcs)
    genotypes = TargeneCore.call_genotypes(
        bgen_prefix, 
        Set(string.(variables.genetic_variants)), 
        call_threshold
    )
    data = TargeneCore.merge(traits, pcs, genotypes)

    # Parameter files
    estimands = TargeneCore.adjusted_parameters(
        estimands, variables, data; 
        positivity_constraint=positivity_constraint
    )

    # write data and parameter files
    write_tmle_inputs(outprefix, data, estimands, batch_size=batch_size)
end