actors_error() = throw(ArgumentError("At least two of (bqtls, env, trans-actors) should be specified"))

treatments_from_actors(::Nothing, ::Nothing, _) = actors_error()
treatments_from_actors(_, ::Nothing, ::Nothing) = actors_error()
treatments_from_actors(::Nothing, _, ::Nothing) = actors_error()

function treatments_from_actors(bqtl_file, env_file, trans_actors_prefix)
    bqtls = read_snps_from_csv(bqtl_file)
    extraT  = read_txt_file(env_file)
    transactors = trans_actors_from_prefix(trans_actors_prefix)
    bqtls, transactors, extraT
end


combine_trans_actors(trans_actors::Vector{DataFrame}, extraT::DataFrame, order) = combinations([trans_actors..., extraT], order)
combine_trans_actors(trans_actors::Vector{DataFrame}, extraT::Nothing, order) = combinations(trans_actors, order)
combine_trans_actors(trans_actors::Nothing, extraT::DataFrame, order) = [[extraT]]


function combine_by_bqtl(bqtls::DataFrame, trans_actors::Union{Vector{DataFrame}, Nothing}, extraT::Union{DataFrame, Nothing}, order::Int)
    param_files = Dict[]
    if order == 1
        for bqtl in bqtls.ID
            push!(param_files, Dict{Any, Any}("T" => [bqtl]))
        end
    else
        for comb in combine_trans_actors(trans_actors, extraT, order - 1)
            comb_interactions = crossjoin(bqtls, comb..., makeunique=true)
            for row in eachrow(comb_interactions)
                push!(
                    param_files,
                    Dict{Any, Any}(
                        "T" => vcat(row["ID"], [row[string("ID_", i)] for i in 1:order-1])
                    )
                )
            end
        end
    end
    return param_files
end

parse_orders(orders::String) = parse.(Int, split(orders, ","))


all_snps(bqtls::DataFrame, transactors::Nothing) = bqtls.ID
all_snps(bqtls::Nothing, transactors::Vector{DataFrame}) = vcat((x.ID for x in transactors)...)
all_snps(bqtls::DataFrame, transactors::Vector{DataFrame}) = vcat(bqtls.ID, (x.ID for x in transactors)...)


read_snps_from_csv(path::Nothing) = nothing
read_snps_from_csv(path::String) = unique(CSV.read(path, DataFrame; select=[:ID, :CHR]), :ID)


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

sorted_unique_values(v) = sort(unique(skipmissing(v)))

"""
    control_case_settings(::Type{IATE}, treatment_tuple, data::DataFrame)

Generates all possible IATE parameter treatment settings where control ≠ case. 
The procedure works as follows:

    - For each treatment, control ≠ case are generated from pairwise combinations of unique values
    - The cartesian product of all those treatment case/control values is returned
"""
function control_case_settings(::Type{IATE}, treatment_tuple, data::DataFrame)
    control_cases_ = []
    # For each treatment variable generate the case/control options
    # where case and control are different
    for t in treatment_tuple
        unique_vals = sorted_unique_values(data[!, t])
        push!(control_cases_, collect(combinations(unique_vals, 2)))
    end
    # Generate the product of all treatment case/control settings
    return Iterators.product(control_cases_...)
end

"""
    control_case_settings(::Type{ATE}, treatment_tuple, data)

Generates all possible ATE parameter treatment settings where only the bQTL has control ≠ case, 
for the rest of the treatment variables: case=control. The bQTl is assumed to 
be the first variable in the `treatment_tuple`.

The procedure works as follows:

    - For the bQTL, control ≠ case are generated from pairwise combinations of unique alleles
    - For all other treatment variables, control = case are the unique values
    - The cartesian product of all those treatment case/control values is returned
"""
function control_case_settings(::Type{ATE}, treatment_tuple, data)
    control_cases_ = []
    # The bQTL is the first treatment variable
    # We get all case/control values
    bQTL_values = TargeneCore.sorted_unique_values(data[!, treatment_tuple[1]])
    push!(control_cases_, collect(combinations(bQTL_values, 2)))
    # The remainder of the treatment are kep fixed
    for t in treatment_tuple[2:end]
        unique_vals = TargeneCore.sorted_unique_values(data[!, t])
        push!(control_cases_, [[uv, uv] for uv in unique_vals])
    end
    # Generate the product of all treatment case/control settings
    return Iterators.product(control_cases_...)
end

function addParameter!(params, param_type::Type{<:TMLE.Parameter}, setting, treatment_tuple, freqs; positivity_constraint=0.)
    if TargeneCore.satisfies_positivity(param_type, setting, freqs; positivity_constraint=positivity_constraint)
        param = Dict{String, Any}("name" => replace(string(param_type), "TMLE." => ""))
        for var_index in eachindex(treatment_tuple)
            control, case = setting[var_index]
            treatment = treatment_tuple[var_index]
            param[treatment] = Dict("control" => control, "case" => case)
        end
        push!(params, param)
    end
end

function addParameters!(param_file, data; positivity_constraint=0.)
    treatment_tuple = param_file["T"]
    freqs = TargeneCore.frequency_table(data, treatment_tuple)
    param_file["Parameters"] = Dict[]
    # This loop adds all ATE parameters where all other treatments than
    # the bQTL are fixed, at the order 1, this is the simple bQTL's ATE
    for setting in control_case_settings(ATE, treatment_tuple, data)
        addParameter!(param_file["Parameters"], ATE, setting, treatment_tuple, freqs; positivity_constraint=positivity_constraint)
    end
    # This loop adds all IATE parameters that pass the positivity threshold
    if size(treatment_tuple, 1) >= 2
        for setting in control_case_settings(IATE, treatment_tuple, data)
            addParameter!(param_file["Parameters"], IATE, setting, treatment_tuple, freqs; positivity_constraint=positivity_constraint)
        end
    end
end

function get_variables(pcs, traits, extraW, extraC, extraT)
    W = all_confounders(pcs, extraW)
    nontargets = Set(vcat("SAMPLE_ID", extraW, extraC, extraT))
    return Dict(
        "PCs" => pcnames(pcs),
        "extraT" => extraT,
        "extraC" => extraC,
        "extraW" => extraW,
        "W" => W,
        "Y" => targets_from_traits(traits, nontargets)
    )
end

function param_files_from_actors(bqtls, transactors, data, variables, orders; batch_size=nothing, positivity_constraint=0.)
    new_param_files = Dict[]
    # For each interaction order generate parameter files
    for order in orders
        # First generate the `T` section
        extraT_df = variables["extraT"] isa Nothing ? nothing : DataFrame(ID=variables["extraT"])
        param_files = TargeneCore.combine_by_bqtl(bqtls, transactors, extraT_df, order)
        for param_file in param_files
            # Generate `Parameters` section
            addParameters!(param_file, data; positivity_constraint=positivity_constraint)
            # If at least one parameter satisfies positivity
            if param_file["Parameters"] != []
                # Generate `W` section
                param_file["W"] = variables["W"]
                # Generate `C` section
                if !(variables["extraC"] isa Nothing)
                    param_file["C"] = variables["extraC"]
                end
                
                add_batchified_param_files!(new_param_files, param_file, variables["Y"], batch_size)
            end
        end
    end
    return new_param_files
end


function tmle_inputs_from_actors(parsed_args)
    batch_size = parsed_args["phenotype-batch-size"]
    outprefix = parsed_args["out-prefix"]
    call_threshold = parsed_args["call-threshold"]
    bgen_prefix = parsed_args["bgen-prefix"]
    positivity_constraint = parsed_args["positivity-constraint"]
    traits = TargeneCore.read_data(parsed_args["traits"])
    pcs = TargeneCore.read_data(parsed_args["pcs"])
    orders = TargeneCore.parse_orders(parsed_args["from-actors"]["orders"])
    extraW = TargeneCore.read_txt_file(parsed_args["from-actors"]["extra-confounders"])
    extraC = TargeneCore.read_txt_file(parsed_args["from-actors"]["extra-covariates"])
    genotypes_asint = parsed_args["from-actors"]["genotypes-as-int"]

    # Retrieve SNPs and environmental treatments
    bqtls, transactors, extraT = TargeneCore.treatments_from_actors(
        parsed_args["from-actors"]["bqtls"], 
        parsed_args["from-actors"]["extra-treatments"], 
        parsed_args["from-actors"]["trans-actors-prefix"]
    )
    # Genotypes and final dataset
    snp_list = TargeneCore.all_snps(bqtls, transactors)
    genotypes = TargeneCore.call_genotypes(bgen_prefix, snp_list, call_threshold; asint=genotypes_asint)
    data = TargeneCore.merge(traits, pcs, genotypes)

    # Parameter files
    variables = TargeneCore.get_variables(pcs, traits, extraW, extraC, extraT)
    param_files = TargeneCore.param_files_from_actors(
        bqtls, transactors, data, variables, orders; 
        batch_size=batch_size, positivity_constraint=positivity_constraint
    )

    # write data and parameter files
    write_tmle_inputs(outprefix, data, param_files)
end