
"""
    instantiate_dataset(path::String)

Returns a DataFrame wrapper around a dataset, either in CSV format.
"""
instantiate_dataset(path::String) =
    endswith(path, ".csv") ? CSV.read(path, DataFrame, ntasks=1) : DataFrame(Arrow.Table(path))

BGEN.rsid(s::Symbol) = s

treatment_variables(Ψ::ComposedEstimand) =
    unique(vcat((treatment_variables(arg) for arg ∈ Ψ.args)...))

treatment_variables(Ψ) = collect(keys(Ψ.treatment_values))

outcome_variables(Ψ) = [Ψ.outcome]

outcome_variables(Ψ::ComposedEstimand) = 
    unique(vcat((outcome_variables(arg) for arg ∈ Ψ.args)...))

"""
    getconfounders(v)

Split string and remove principal components from the list.
"""
getconfounders(v) = Symbol.(filter(x -> !occursin(r"^PC[0-9]*$", x), split_string(v)))

is_significant(Ψ̂::TMLE.Estimate; threshold=0.05) = 
    pvalue(OneSampleTTest(Ψ̂)) < threshold

function is_significant(Ψ̂::TMLE.ComposedEstimate; threshold=0.05)
    sig = if length(Ψ̂.estimate) > 1 
        pvalue(TMLE.OneSampleHotellingT2Test(Ψ̂)) < threshold
    else 
        pvalue(TMLE.OneSampleTTest(Ψ̂)) < threshold
    end
    return sig
end

"""
For FailedEstimates
"""
is_significant(Ψ̂; threshold=0.05) = false

"""
For NamedTuples/Dicts stored in a results file
"""
is_significant(nt, estimator_key; threshold=0.05) = 
    is_significant(nt[estimator_key]; threshold=threshold)

function read_significant_from_hdf5(filename; threshold=0.05, estimator_key=:TMLE)
    jldopen(filename) do io
        return mapreduce(vcat, keys(io)) do key
            [nt[estimator_key].estimand for nt ∈ io[key] if is_significant(nt, estimator_key, threshold=threshold)]
        end
    end
end

function read_significant_from_jls(filename; threshold=0.05, estimator_key=:TMLE)
    results = []
    open(filename) do io
        while !eof(io)
            nt = deserialize(io)
            if is_significant(nt, estimator_key; threshold=threshold)
                push!(results, nt[estimator_key].estimand)
            end
        end
    end
    return results
end

read_significant_from_json(filename; threshold=0.05, estimator_key=:TMLE) = 
    [nt[estimator_key].estimand for nt ∈ TMLE.read_json(filename) if is_significant(nt, estimator_key; threshold=threshold)]

function read_significant_results(filename; threshold=0.05, estimator_key=:TMLE)
    results = if endswith(filename, "hdf5")
        read_significant_from_hdf5(filename; threshold=threshold, estimator_key=estimator_key)
    elseif endswith(filename, "jls")
        read_significant_from_jls(filename; threshold=threshold, estimator_key=estimator_key)
    elseif endswith(filename, "json")
        read_significant_from_json(filename; threshold=threshold, estimator_key=estimator_key)
    else
        throw(ArgumentError("Unupported estimate file format: $filepath"))
    end
    return results
end

"""
    group_parameters(parameters; min_size=1)

Tries to group parameters together in optimal ordering 
by group of size approximately greater than min_size.
"""
function write_estimands_files(prefix, parameters, chunksize; fileprefix="estimands_")
    for (batch_index, batch) in enumerate(Iterators.partition(parameters, chunksize))
        filename = string(fileprefix, batch_index, ".jls")
        serialize(filepath_from_prefix(prefix; filename=filename), Configuration(estimands=batch))
    end
end

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

write_dataset(prefix, dataset; filename="dataset.arrow") = 
    Arrow.write(filepath_from_prefix(prefix; filename=filename), dataset)

set_from_txt_file(filepath::AbstractString) = Set(open(readlines, filepath))
