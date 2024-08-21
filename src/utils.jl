###############################################################################
###                             MISC UTILITIES                              ###
###############################################################################

"""
    read_csv_file(filepath)

The SAMPLE_ID column should be read as a String.
"""
read_csv_file(filepath) = CSV.read(filepath, DataFrame, types=Dict(:SAMPLE_ID => String))

function files_matching_prefix(prefix)
    directory, _prefix = splitdir(prefix)
    _directory = directory == "" ? "." : directory

    return map(
        f -> joinpath(directory, f),
        filter(
            f -> startswith(f, _prefix), 
            readdir(_directory)
        )
    )
end

function make_filepath_from_prefix(prefix; filename="dataset.arrow")
    return if isdir(prefix)
        joinpath(prefix, filename)
    else
        dir, _prefix = splitdir(prefix)
        joinpath(dir, string(_prefix, filename))
    end
end

function merge(traits, pcs, genotypes)
    return innerjoin(
        innerjoin(traits, pcs, on="SAMPLE_ID"),
        genotypes,
        on="SAMPLE_ID"
    )
end

###############################################################################
###                      ESTIMATION INPUTS WRITING                          ###
###############################################################################

batch_name(outprefix, batch_id) = string(outprefix, ".estimands_", batch_id, ".jls") 

save_estimands(outprefix, estimands, batchsize::Nothing) = 
    serialize(batch_name(outprefix, 1), Configuration(estimands=estimands))

function save_estimands(outprefix, estimands, batchsize)
    for (batch_id, batch) ∈ enumerate(Iterators.partition(estimands, batchsize))
        batchfilename = batch_name(outprefix, batch_id)
        serialize(batchfilename, Configuration(estimands=batch))
    end
end

function write_estimation_inputs(outprefix, final_dataset, estimands; batchsize=nothing)
    # Write final_dataset
    Arrow.write(string(outprefix, ".data.arrow"), final_dataset)
    # Write estimands_files
    save_estimands(outprefix, estimands, batchsize)
end

###############################################################################
###                      GENOTYPES CALLING FROM BGEN                        ###
###############################################################################

const CHR_REG = r"chr[1-9]+"

NotAllVariantsFoundError(rsids) = 
    ArgumentError(string("Some variants were not found in the genotype files: ", join(rsids, ", ")))

NotBiAllelicOrUnphasedVariantError(rsid) = ArgumentError(string("Variant: ", rsid, " is not bi-allelic or not unphased."))

BGEN.rsid(s::Symbol) = s

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

"""
    genotypes_encoding(variant)

String genotypes are reported.
"""
function genotypes_encoding(variant)
    all₁, all₂ = alleles(variant)
    return [all₁*all₁, all₁*all₂, all₂*all₂]
end

call_genotype(variant_genotypes, max_index, max_prob, threshold::Nothing) = variant_genotypes[max_index]

call_genotype(variant_genotypes, max_index, max_prob, threshold::Real) = max_prob >= threshold ? variant_genotypes[max_index] : missing

function call_genotypes(probabilities::AbstractArray, variant_genotypes::AbstractVector{T}, threshold) where T
    n = size(probabilities, 2)
    t = Vector{Union{T, Missing}}(missing, n)
    for i in 1:n
        # If no allele has been annotated with sufficient confidence
        # the sample is declared as missing for this variant
        max_prob, max_index = findmax(probabilities[:, i])
        if !isinf(max_prob)
            t[i] = call_genotype(variant_genotypes, max_index, max_prob, threshold)
        end
    end
    return t
end

"""
    bgen_files(snps, bgen_prefix)

This function assumes the UK-Biobank structure
"""
function call_genotypes(bgen_prefix::String, query_rsids::Set{<:AbstractString}, threshold)
    query_rsids = copy(query_rsids)
    chr_dir_, prefix_ = splitdir(bgen_prefix)
    chr_dir = chr_dir_ == "" ? "." : chr_dir_
    genotypes = nothing
    for filename in readdir(chr_dir)
        length(query_rsids) == 0 ? break : nothing
        if is_numbered_chromosome_file(filename, prefix_)
            bgenfile = read_bgen(joinpath(chr_dir_, filename))
            chr_genotypes = DataFrame(SAMPLE_ID=bgenfile.samples)
            bgen_rsids = Set(rsids(bgenfile))
            for query_rsid in query_rsids
                if query_rsid ∈ bgen_rsids
                    pop!(query_rsids, query_rsid)
                    variant = variant_by_rsid(bgenfile, query_rsid)
                    if n_alleles(variant) != 2
                        @warn("Skipping $query_rsid, not bi-allelic")
                        continue
                    end
                    minor_allele_dosage!(bgenfile, variant)
                    variant_genotypes = genotypes_encoding(variant)
                    probabilities = probabilities!(bgenfile, variant)
                    probabilities[isnan.(probabilities)] .= -Inf
                    size(probabilities, 1) != 3 && throw(NotBiAllelicOrUnphasedVariantError(query_rsid))
                    chr_genotypes[!, query_rsid] = call_genotypes(probabilities, variant_genotypes, threshold)
                end
            end
            genotypes = genotypes isa Nothing ? chr_genotypes :
                    innerjoin(genotypes, chr_genotypes, on=:SAMPLE_ID)
        end
    end
    length(query_rsids) == 0 || throw(NotAllVariantsFoundError(query_rsids))
    return genotypes
end

###############################################################################
###                          ESTIMANDS ACCESSORS                            ###
###############################################################################

MismatchedVariableError(variable) = ArgumentError(string("Each component of a JointEstimand should contain the same ", variable, " variables."))

function check_variables_are_consistent(Ψ::JointEstimand, variable, variable_accessor, args...)
    if length(Ψ.args) > 1
        for Ψᵢ in Ψ.args[2:end]
            variable_accessor(Ψᵢ, args...) == variable || throw(MismatchedVariableError(split(string(variable_accessor), "_")[end]))
        end
    end
end

get_outcome(Ψ) = Ψ.outcome

function get_outcome(Ψ::JointEstimand)
    outcome = get_outcome(first(Ψ.args))
    check_variables_are_consistent(Ψ, outcome, get_outcome)
    return outcome
end

# Updated for OrderedDict
get_treatments(Ψ) = sort(collect(keys(Ψ.treatment_values)))

function get_treatments(Ψ::JointEstimand)
    treatments = get_treatments(first(Ψ.args))
    check_variables_are_consistent(Ψ, treatments, get_treatments)
    return treatments
end

function get_all_confounders(Ψ::JointEstimand)
    all_confounders = get_all_confounders(first(Ψ.args))
    check_variables_are_consistent(Ψ, all_confounders, get_all_confounders)
    return all_confounders
end

# Updated for OrderedDict
get_all_confounders(Ψ) = Tuple(sort(collect(Iterators.flatten(( value for (key, value) ∈ pairs(Ψ.treatment_confounders))))))

get_confounders(Ψ, treatment) = Ψ.treatment_confounders[treatment]

function get_confounders(Ψ::JointEstimand, treatment)
    confounders = get_confounders(first(Ψ.args), treatment)
    check_variables_are_consistent(Ψ, confounders, get_confounders, treatment)
    return confounders
end

get_outcome_extra_covariates(Ψ) = Ψ.outcome_extra_covariates

function get_outcome_extra_covariates(Ψ::JointEstimand)
    outcome_extra_covariates = get_outcome_extra_covariates(first(Ψ.args))
    check_variables_are_consistent(Ψ, outcome_extra_covariates, get_outcome_extra_covariates)
    return outcome_extra_covariates
end

###############################################################################
###                    P-VALUES WITH ERROR MANAGEMENT                       ###
###############################################################################

default_null(Ψ̂::TMLE.JointEstimate) = zeros(length(Ψ̂.estimates))

default_null(Ψ̂::TMLE.EICEstimate) = 0.

function pvalue_or_nan(Ψ̂, Ψ₀=default_null(Ψ̂))
    return try
        pvalue(significance_test(Ψ̂, Ψ₀))
    catch
        NaN
    end
end

pvalue_or_nan(Ψ̂::TmleCLI.FailedEstimate, Ψ₀=nothing) = NaN