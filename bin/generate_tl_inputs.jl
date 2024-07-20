using TargeneCore
using ArgParse

function parse_commandline()
    s = ArgParseSettings(description="Preparation of data for TMLE, two modes are currently available.")

    @add_arg_table s begin
        "from-param-file"
            action = :command
            help = "Will assume Parameter files (see README.md) to be given."

        "allele-independent"
            action = :command
            help = "Generates allele independent estimands from a configuration file."

        "--out-prefix"
            arg_type = String
            help = "Output prefix"
            default = "final"

        "--batch-size"
            arg_type = Int
            required = false
            help = "To split estimands in multiple files of size batch-size"

        "--positivity-constraint"
            arg_type = Float64
            default = 0.01
            required = false
            help = "Minimum frequency a treatment value occuring in a Parameter must reach in the population"

        "--verbosity", "-v"
            help = "Verbosity level"
            arg_type = Int
            default = 1
    end

    @add_arg_table s["from-param-file"] begin
        "paramfile"
            arg_type = String
            help = "Parameter file"
        
        "--traits"
            arg_type = String
            required = true
            help = "Path to the traits dataset"

        "--call-threshold"
            arg_type = Float64
            help = "Hard call threshold for imputed genotypes"
            default = nothing

        "--bgen-prefix"
            help = "Prefix path to BGEN chromosomes."
            required = true
            arg_type = String

        "--pcs"
            arg_type = String
            required = true
            help = "Path to a genetic confounders file"
    end

    @add_arg_table s["allele-independent"] begin
        "config"
            arg_type = String
            help = "Configuration file (.yaml)."
        
        "--traits"
            arg_type = String
            required = true
            help = "Path to the traits dataset"

        "--call-threshold"
            arg_type = Float64
            help = "Hard call threshold for imputed genotypes"
            default = nothing

        "--genotype-prefix"
            help = "Prefix path to BGEN/BED chromosomes."
            required = true
            arg_type = String

        "--pcs"
            arg_type = String
            required = true
            help = "Path to a genetic confounders file"
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

tl_inputs(parsed_args)