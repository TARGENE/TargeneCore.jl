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

BGEN.rsid(s::Symbol) = s

treatment_variables(Ψ::JointEstimand) =
    unique(vcat((treatment_variables(arg) for arg ∈ Ψ.args)...))

treatment_variables(Ψ) = collect(keys(Ψ.treatment_values))

outcome_variables(Ψ) = [Ψ.outcome]

outcome_variables(Ψ::JointEstimand) = 
    unique(vcat((outcome_variables(arg) for arg ∈ Ψ.args)...))

"""
    getconfounders(v)

Split string and remove principal components from the list.
"""
getconfounders(v) = Symbol.(filter(x -> !occursin(r"^PC[0-9]*$", x), split_string(v)))

function unique_treatments(estimands)
    treatments = Set{String}()
    for Ψ in estimands
        union!(treatments, string.(treatment_variables(Ψ)))
    end
    return treatments
end

function filepath_from_prefix(prefix; filename="dataset.arrow")
    return if isdir(prefix)
        joinpath(prefix, filename)
    else
        dir, _prefix = splitdir(prefix)
        joinpath(dir, string(_prefix, filename))
    end
end

set_from_txt_file(filepath::AbstractString) = Set(open(readlines, filepath))

#### From tl_inputs

const CHR_REG = r"chr[1-9]+"

batch_name(outprefix, batch_id) = string(outprefix, ".estimands_", batch_id, ".jls") 

"""
    read_csv_file(filepath)

The SAMPLE_ID column should be read as a String.
"""
read_csv_file(filepath) = CSV.read(filepath, DataFrame, types=Dict(:SAMPLE_ID => String))

function write_tl_inputs(outprefix, final_dataset, estimands; batch_size=nothing)
    # Write final_dataset
    Arrow.write(string(outprefix, ".data.arrow"), final_dataset)
    # Write param_files
    if batch_size !== nothing
        for (batch_id, batch) in enumerate(Iterators.partition(estimands, batch_size))
            serialize(batch_name(outprefix, batch_id), Configuration(estimands=batch))
        end
    else
        serialize(string(outprefix, ".estimands.jls"), Configuration(estimands=estimands))
    end
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

NotAllVariantsFoundError(rsids) = 
    ArgumentError(string("Some variants were not found in the genotype files: ", join(rsids, ", ")))

NotBiAllelicOrUnphasedVariantError(rsid) = ArgumentError(string("Variant: ", rsid, " is not bi-allelic or not unphased."))
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

read_txt_file(path::Nothing) = nothing
read_txt_file(path) = CSV.read(path, DataFrame, header=false)[!, 1]

pcnames(pcs) = filter(!=("SAMPLE_ID"), names(pcs))

all_confounders(pcs, extraW::Nothing) = Symbol.(pcnames(pcs))
all_confounders(pcs, extraW) = Symbol.(unique(vcat(pcnames(pcs), extraW)))

targets_from_traits(traits, non_targets) = Symbol.(filter(x -> x ∉ non_targets, names(traits)))


function merge(traits, pcs, genotypes)
    return innerjoin(
        innerjoin(traits, pcs, on="SAMPLE_ID"),
        genotypes,
        on="SAMPLE_ID"
    )
end

estimand_with_new_outcome(Ψ::T, outcome) where T = T(
    outcome=outcome, 
    treatment_values=Ψ.treatment_values, 
    treatment_confounders=Ψ.treatment_confounders, 
    outcome_extra_covariates=Ψ.outcome_extra_covariates
)

function update_estimands_from_outcomes!(estimands, Ψ::T, outcomes) where T
    for outcome in outcomes
        push!(
            estimands, 
            T(
                outcome=outcome, 
                treatment_values=Ψ.treatment_values, 
                treatment_confounders=Ψ.treatment_confounders, 
                outcome_extra_covariates=Ψ.outcome_extra_covariates)
        )
    end
end


function update_estimands_from_outcomes!(estimands, Ψ::JointEstimand, outcomes)
    for outcome in outcomes
        push!(
            estimands,
            JointEstimand((estimand_with_new_outcome(arg, outcome) for arg in Ψ.args)...)
        )
    end
end
