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

function call_genotypes(probabilities::AbstractArray, variant_genotypes::AbstractVector{T}, threshold::Real) where T
    n = size(probabilities, 2)
    t = Vector{Union{T, Missing}}(missing, n)
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
function call_genotypes(bgen_prefix::String, query_rsids::Set{<:AbstractString}, threshold::Real)
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

sorted_treatment_names(Ψ) = tuple(sort(collect(keys(Ψ.treatment_values)))...)

function setting_iterator(Ψ::TMLE.StatisticalIATE)
    treatments = sorted_treatment_names(Ψ)
    return (
        NamedTuple{treatments}(collect(Tval)) for 
            Tval in Iterators.product((values(Ψ.treatment_values[T]) for T in treatments)...)
    )
end

function setting_iterator(Ψ::TMLE.StatisticalATE)
    treatments = sorted_treatment_names(Ψ)
    return (
        NamedTuple{treatments}([(Ψ.treatment_values[T][c]) for T in treatments])
            for c in (:case, :control)
    )
end

function setting_iterator(Ψ::TMLE.StatisticalCM)
    treatments = sorted_treatment_names(Ψ)
    return (NamedTuple{treatments}(Ψ.treatment_values[T] for T in treatments), )
end

function satisfies_positivity(Ψ::TMLE.Estimand, freqs; positivity_constraint=0.01)
    for base_setting in setting_iterator(Ψ)
        if !haskey(freqs, base_setting) || freqs[base_setting] < positivity_constraint
            return false
        end
    end
    return true
end

function frequency_table(data, treatments::AbstractVector)
    treatments = sort(treatments)
    freqs = Dict()
    N = nrow(data)
    for (key, group) in pairs(groupby(data, treatments; skipmissing=true))
        freqs[NamedTuple(key)] = nrow(group) / N
    end
    return freqs
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

"""
    tl_inputs(parsed_args)

Support for the generation of estimands according to 2 strategies:
- from-actors
- from-param-file
"""
function tl_inputs(parsed_args)
    if parsed_args["%COMMAND%"] == "from-actors"
        tl_inputs_from_actors(parsed_args)
    elseif parsed_args["%COMMAND%"] == "from-param-file"
        tl_inputs_from_param_files(parsed_args)
    elseif parsed_args["%COMMAND%"] == "allele-independent"
        allele_independent_estimands(parsed_args)
    elseif parsed_args["%COMMAND%"] == "permutation-tests"
        permutation_tests_tl_inputs(parsed_args)
    elseif parsed_args["%COMMAND%"] == "loco-gwas"
        loco_gwas(parsed_args)
    else
        throw(ArgumentError("Unrecognized command."))
    end
end


