const ENSEMBL_URL = "https://rest.ensembl.org"
const REGULATORY_ELEMENTS = Set(["exon"])

NotEnoughMatchingVariantsError(rsid, p, reltol) = 
    ArgumentError(string("Could not find ", p, 
    " random matching variants for ", rsid, ", consider decreasing the MAF reltol (current: ", reltol, ")"
    ))

"""
    update_treatment_setting!(variant::Variant, origin_variant::Variant, origin_case, origin_control)

When the origin_variant is a `Variant` it means it is a "mapped" transactor and
new alleles must be inferred.
"""
function get_new_treatment_setting(variant::Variant, origin_variant::Variant, origin_setting)
    allele_map = Dict(
        minor_allele(origin_variant) => minor_allele(variant),
        major_allele(origin_variant) => major_allele(variant)
    )
    control = join(allele_map[string(a)] for a in origin_setting.control)
    case = join(allele_map[string(a)] for a in origin_setting.case)
    return (control=control, case=case)
end

"""
    update_treatment_setting!(treatment_setting, treatment, origin_treatment, origin_case, origin_control)

When the origin_treatment is anything else it means it is a "non-mapped" treatment and treatment settings are kept identical.
"""
get_new_treatment_setting(treatment, origin_treatment, origin_setting) = origin_setting

"""
    overlaps(v::Variant)

Checks whether a variant is in a regulatory region as defined by ENSEMBL. nm,./\
"""
function is_in_regulatory_region(v::Variant)
    chr = parse(Int, chrom(v))
    loc = pos(v)
    ext = string("/overlap/region/human/", chr, ":", loc, "-", loc, "?feature=regulatory;feature=exon;feature=other_regulatory")
    url = string(ENSEMBL_URL, ext)
    resp = HTTP.get(url, headers=Dict( "Content-Type" => "application/json"))
    if length(JSON.parse(String(resp.body))) == 0
        return false
    else
        return true
    end
end

function same_maf(b::Bgen, variant::Variant, maf::AbstractFloat; reltol=0.05)
    variant_maf = mean(minor_allele_dosage!(b, variant))
    return abs(maf - variant_maf)/maf <= reltol
end

isSNP(v::Variant) = all(length(a) == 1 for a in alleles(v))

function matches_constraints!(b::Bgen, candidate_variant::Variant, variant::Variant, target_maf::AbstractFloat, rs_ids::Set{<:AbstractString}; reltol=0.05)
    rsid(candidate_variant) ∈ rs_ids && return false
    if isSNP(variant)
        isSNP(candidate_variant) || return false
    end
    if same_maf(b, candidate_variant, target_maf; reltol=reltol)
        isreg = false
        try
            isreg = is_in_regulatory_region(candidate_variant)
        catch
            return false
        end
        return !isreg
    end
    return false
end

"""
find_maf_matching_random_variants!(b::Bgen, variant::Variant, rs_ids::Set{<:AbstractString}; p=10, rng=StableRNG(123), reltol=0.05)

For the given variant find `p` random variants in the BGEN file
that match the variant's MAF (up to the relative `tol`) and are not in:

- (i) exons
- (ii) promoters
- (iii) enhancers
- (iv) `rs_ids`

"""
function find_maf_matching_random_variants!(
    b::Bgen, variant::Variant, rs_ids::Set{<:AbstractString}; 
    p=10, rng=StableRNG(123), reltol=0.05, verbosity=1
    )
    target_maf = mean(minor_allele_dosage!(b, variant))
    all_rsids = shuffle(rng, rsids(b))
    matching_variants = Vector{Variant}(undef, p)
    index = 1
    for candidate_rs_id in all_rsids
        # Some rs_ids do not correspond to a unique entry and can raise,
        # we just skip those
        candidate_variant = try
            variant_by_rsid(b, candidate_rs_id)
        catch
            verbosity > 0 && @info(string(
                "An error occured while looking for variant ", 
                candidate_rs_id,
                " in BGEN file: continuing."))
            continue
        end
        if matches_constraints!(b, candidate_variant, variant, target_maf, rs_ids; reltol=reltol)
            verbosity > 0 && @info("Found matching variant: ", rsid(candidate_variant))
            matching_variants[index] = candidate_variant
            index += 1
            index > p && return matching_variants
        end
    end
    throw(NotEnoughMatchingVariantsError(rsid(variant), p, reltol))
end

function find_maf_matching_random_variants(
    trans_actors::Set{<:AbstractString}, bgen_prefix::AbstractString, all_rsids::Set{<:AbstractString}; 
    p=10, rng=StableRNG(123), reltol=0.05, verbosity=1
    )
    chr_dir_, prefix_ = splitdir(bgen_prefix)
    chr_dir = chr_dir_ == "" ? "." : chr_dir_
    variant_map = Dict()
    remaining_transactors = deepcopy(trans_actors)
    for filename in readdir(chr_dir)
        length(remaining_transactors) == 0 && return variant_map
        if is_numbered_chromosome_file(filename, prefix_)
            b = read_bgen(joinpath(chr_dir_, filename))
            bgen_rsids = Set(rsids(b))
            for rs_id in remaining_transactors
                if rs_id ∈ bgen_rsids
                    variant = variant_by_rsid(b, rs_id)
                    verbosity > 0 && @info("Looking for variants matching: ", rs_id)
                    minor_allele_dosage!(b, variant)
                    variant_map[Symbol(rs_id)] = (
                        variant,
                        find_maf_matching_random_variants!(
                        b, variant, all_rsids;
                        p=p, rng=rng, reltol=reltol, verbosity=verbosity)
                    )
                    pop!(remaining_transactors, rs_id)
                end
            end
        end
    end
    return variant_map
end


function make_random_variants_estimands(Ψ::ComposedEstimand, variant_map)
    newargs = Tuple(make_random_variants_estimands(arg, variant_map) for arg ∈ Ψ.args)
    return [ComposedEstimand(Ψ.f, Tuple(args)) for args ∈ zip(newargs...)]
end

function make_random_variants_estimands(Ψ::T, variant_map) where T <: TMLE.Estimand
    transactors = keys(variant_map)
    origin_treatment_variables = keys(Ψ.treatment_values)
    # At least one trans-actor in the parameter treatments to be processed
    new_estimands = []
    if any(rs_id ∈ origin_treatment_variables for rs_id in transactors)
        new_treatments = [t ∈ transactors ? variant_map[t][2] : [t] for t in origin_treatment_variables]
        for new_treatment_vars ∈ Iterators.product(new_treatments...)
            new_treatment_names = Tuple(Symbol(rsid(v)) for v in new_treatment_vars)
            new_treatments_setting = []
            for (origin_treatment_name, new_treatment_var) in zip(origin_treatment_variables, new_treatment_vars)
                # Potentially retrieve the genetic variant to match alleles
                origin_treatment_var = haskey(variant_map, origin_treatment_name) ? variant_map[origin_treatment_name][1] : new_treatment_var
                # Update treatment settings with new matched alleles
                new_treatment_setting = get_new_treatment_setting(
                    new_treatment_var, 
                    origin_treatment_var, 
                    Ψ.treatment_values[origin_treatment_name]
                )
                push!(new_treatments_setting, new_treatment_setting)
            end
            treatment_values = NamedTuple{new_treatment_names}(new_treatments_setting)
            treatment_confounders = NamedTuple{new_treatment_names}(values(Ψ.treatment_confounders))
            
            push!(
                new_estimands,
                T(;
                    outcome = Ψ.outcome,
                    treatment_values = treatment_values,
                    treatment_confounders = treatment_confounders,
                    outcome_extra_covariates = Ψ.outcome_extra_covariates
                )
            )
        end
    end
    return new_estimands
end

make_random_variants_estimands(estimands, variant_map) = 
    vcat((make_random_variants_estimands(Ψ, variant_map) for Ψ in estimands)...)

function generate_random_variants_estimands(parsed_args)
    resultsfile = parsed_args["results"]
    p = parsed_args["p"]
    rng = StableRNG(parsed_args["rng"])
    reltol = parsed_args["reltol"]
    trans_actors = set_from_txt_file(parsed_args["variants-to-randomize"])
    bgen_prefix = parsed_args["bgen-prefix"]
    pval_threshold = parsed_args["pval-threshold"]
    out = parsed_args["out"]
    verbosity = parsed_args["verbosity"]
    estimator_key = Symbol(parsed_args["estimator-key"])
    
    verbosity > 0 && @info string("Retrieving significant estimands.")
    significant_estimands = read_significant_results(resultsfile; threshold=pval_threshold, estimator_key=estimator_key)
    all_rsids = unique_treatments(significant_estimands)

    verbosity > 0 && @info string("Looking for random MAF matching variants for each Trans-Actor.")
    variant_map = find_maf_matching_random_variants(
        trans_actors, bgen_prefix, all_rsids; 
        p=p, rng=rng, reltol=reltol, verbosity=verbosity
    )
    verbosity > 0 && @info string("Building new estimands from matched random variants.")
    new_estimands = make_random_variants_estimands(significant_estimands, variant_map)
    new_estimands = groups_ordering(new_estimands, brute_force=false, do_shuffle=false)
    serialize(out, Configuration(estimands = new_estimands))
    verbosity > 0 && @info string("Done.")
    return 0
end