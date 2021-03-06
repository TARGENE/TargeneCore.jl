

using TMLEEpistasis
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Compute the Sieve Variance Plateau estimate for each phenotype in the result file")

    @add_arg_table s begin
        "prefix"
            help = "Prefix to the .hdf5 files generated by TMLEEpistasis `ukbb.jl` script"
            arg_type = String
            required = true
        "grm-prefix"
            arg_type = String
            help = "Prefix of the aggregated GRM"
            required = true
        "out"
            arg_type = String
            help = "output filename"
            required = true
        "--nb-estimators", "-n"
            arg_type = Int
            help = "Number of variance estimators to compute"
            default = 10
        "--max-tau", "-m"
            arg_type = Float64
            help = "Maximum distance of individuals to take into account (maximum=2)"*
                   "It was witnessed that beyond 0.9, weird limit effects happen"
            default = 0.8
        "--pval", "-p"
            arg_type = Float64
            help = "Only traits/queries pairs that satisfy this pvalue threshold with an independence"*
                    " assumption between individuals will be considered for the sieve."
            default = 0.05
        
    end

    return parse_args(s)
end


parsed_args = parse_commandline()

sieve_variance_plateau(parsed_args)
