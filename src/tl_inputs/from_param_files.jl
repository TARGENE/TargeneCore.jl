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

MismatchedVariableError(variable) = ArgumentError(string("Each component of a ComposedEstimand should contain the same ", variable, " variables."))

function check_genotypes_encoding(val::NamedTuple, type)
    if !(typeof(val.case) <: type && typeof(val.control) <: type)
        throw(MismatchedCaseControlEncodingError())
    end
end

check_genotypes_encoding(val::T, type) where T = 
    T <: type || throw(MismatchedCaseControlEncodingError())


get_treatments(Ψ) = keys(Ψ.treatment_values)

function get_treatments(Ψ::JointEstimand)
    treatments = get_treatments(first(Ψ.args))
    if length(Ψ.args) > 1
        for arg in Ψ.args[2:end]
            get_treatments(arg) == treatments || throw(MismatchedVariableError("treatments"))
        end
    end
    return treatments
end


function get_all_confounders(Ψ::JointEstimand)
    confounders = get_all_confounders(first(Ψ.args))
    if length(Ψ.args) > 1
        for arg in Ψ.args[2:end]
            get_all_confounders(arg) == confounders || throw(MismatchedVariableError("confounders"))
        end
    end
    return confounders
end

get_all_confounders(Ψ) = Tuple(sort(collect(Iterators.flatten((Tconf for Tconf ∈ Ψ.treatment_confounders)))))

get_confounders(Ψ, treatment) = Ψ.treatment_confounders[treatment]

function get_confounders(Ψ::JointEstimand, treatment)
    confounders = get_confounders(first(Ψ.args), treatment)
    if length(Ψ.args) > 1
        for arg in Ψ.args[2:end]
            get_confounders(arg, treatment) == confounders || throw(MismatchedVariableError("confounders"))
        end
    end
    return confounders
end

get_outcome_extra_covariates(Ψ) = Ψ.outcome_extra_covariates

function get_outcome_extra_covariates(Ψ::JointEstimand)
    outcome_extra_covariates = get_outcome_extra_covariates(first(Ψ.args))
    if length(Ψ.args) > 1
        for arg in Ψ.args[2:end]
            get_outcome_extra_covariates(arg) == outcome_extra_covariates || throw(MismatchedVariableError("outcome extra covariates"))
        end
    end
    return outcome_extra_covariates
end

get_outcome(Ψ) = Ψ.outcome

function get_outcome(Ψ::JointEstimand)
    outcome = get_outcome(first(Ψ.args))
    if length(Ψ.args) > 1
        for arg in Ψ.args[2:end]
            get_outcome(arg) == outcome || throw(MismatchedVariableError("outcome"))
        end
    end
    return outcome
end

function get_variables(estimands, traits, pcs)
    genetic_variants = Set{Symbol}()
    others = Set{Symbol}()
    pcs = Set{Symbol}(filter(x -> x != :SAMPLE_ID, propertynames(pcs)))
    alltraits = Set{Symbol}(filter(x -> x != :SAMPLE_ID, propertynames(traits)))
    for Ψ in estimands
        treatments = get_treatments(Ψ)
        confounders = get_all_confounders(Ψ)
        outcome_extra_covariates = get_outcome_extra_covariates(Ψ)
        push!(
            others, 
            outcome_extra_covariates..., 
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

adjust_parameter_sections(Ψ::JointEstimand, variants_alleles, pcs) = 
    JointEstimand((adjust_parameter_sections(arg, variants_alleles, pcs) for arg in Ψ.args)...)

function estimand_satisfying_positivity(Ψ, frequency_table; positivity_constraint=0.01)
    if TMLE.satisfies_positivity(Ψ, frequency_table;positivity_constraint=positivity_constraint)
        return Ψ
    else
        return nothing
    end
end

function estimand_satisfying_positivity(Ψ::JointEstimand, frequency_table; positivity_constraint=0.01)
    new_args = Tuple((arg for arg in Ψ.args if TMLE.satisfies_positivity(arg, frequency_table;positivity_constraint=positivity_constraint)))
    if isempty(new_args)
        return nothing
    else
        return JointEstimand(new_args...)
    end
end
    
function append_from_valid_estimands!(
    estimands::Vector{<:TMLE.Estimand},
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
    treatments = get_treatments(Ψ)
    if !haskey(frequency_tables, treatments)
        frequency_tables[treatments] = TMLE.get_frequency_table(data, treatments)
    end
    # Check estimand or filter joint estimand to satisfy positivity
    Ψ = estimand_satisfying_positivity(Ψ, frequency_tables[treatments]; positivity_constraint=positivity_constraint)
    if !(Ψ === nothing)
        # Expand wildcard to all outcomes
        if get_outcome(Ψ) === :ALL
            update_estimands_from_outcomes!(estimands, Ψ, variables.outcomes)
        else
            # Ψ.target || MissingVariableError(variable)
            push!(estimands, Ψ)
        end
    end
end

function adjusted_estimands(estimands, variables, data; positivity_constraint=0.)
    final_estimands = TMLE.Estimand[]
    variants_alleles = Dict(v => Set(unique(skipmissing(data[!, v]))) for v in variables.genetic_variants)
    freqency_tables = Dict()
    for Ψ in estimands
        # If the genotypes encoding is a string representation make sure they match the actual genotypes
        append_from_valid_estimands!(final_estimands, freqency_tables, Ψ, data, variants_alleles, variables; positivity_constraint=positivity_constraint)
    end

    length(final_estimands) > 0 || throw(NoRemainingParamsError(positivity_constraint))
    ordered_estimands = groups_ordering(final_estimands; brute_force=false, do_shuffle=false)
    
    return ordered_estimands
end


function tl_inputs_from_param_files(config_file, genotypes_prefix, traits_file, pcs_file;
    outprefix="final",
    batchsize=nothing,
    call_threshold=0.9,
    positivity_constraint=0.01,
    verbosity=0)

    # Load initial Parameter files and data
    estimands = TMLE.read_yaml(config_file).estimands
    traits = TargeneCore.read_csv_file(traits_file)
    pcs = TargeneCore.read_csv_file(pcs_file)

    # Genotypes and full data
    verbosity > 0 && @info("Creating dataset.")
    variables = TargeneCore.get_variables(estimands, traits, pcs)
    genotypes = TargeneCore.call_genotypes(
        genotypes_prefix, 
        Set(string.(variables.genetic_variants)), 
        call_threshold
    )
    data = TargeneCore.merge(traits, pcs, genotypes)

    # Adjusted Estimands
    verbosity > 0 && @info("Validating estimands.")
    estimands = TargeneCore.adjusted_estimands(
        estimands, variables, data; 
        positivity_constraint=positivity_constraint
    )

    # write data and parameter files
    write_tl_inputs(outprefix, data, estimands, batch_size=batchsize)
end