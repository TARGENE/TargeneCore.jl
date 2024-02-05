using ArgParse
using TargeneCore

function parse_commandline()
    s = ArgParseSettings(
        description = "Random variant parameters generation",
        commands_are_required = false,
        version = "0.1",
        add_version = true)

    @add_arg_table s begin
        "bgen-prefix"
            help = "Prefix to BGEN files."
            required = true
        "traits"
            help = "File containing traits."
            required = true
        "confounders"
            help = "File containing genetic confounders."
            required = true
        "variants"
            help = "File containing genetic variants to extract."
            required = true
        "--call-threshold"
            arg_type = Float64
            help = "Hard call threshold for imputed genotypes"
            default = 0.9
        "--out"
            help = "Output dataset file."
            arg_type = String
            default = "dataset.arrow"
        "--verbosity"
            help = "Logging level."
            arg_type = Int
            default = 1
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

generate_dataset(parsed_args)

