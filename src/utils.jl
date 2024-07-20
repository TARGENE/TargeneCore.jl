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
