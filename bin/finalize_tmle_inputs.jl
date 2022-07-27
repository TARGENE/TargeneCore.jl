using TMLEEpistasis
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Simply ensures the files have the same order and merges some files.")

    @add_arg_table s begin
        "--out-prefix"
            arg_type = String
            help = "Output prefix"
            default = "final"
        "--genotypes"
            arg_type = String
            help = "Path to a genotypes file"
            required = true
        "--binary-phenotypes"
            arg_type = String
            required = true
            help = help = "Path to the binary phenotypes file"
        "--continuous-phenotypes"
            arg_type = String
            required = true
            help = help = "Path to the continuous phenotypes file"
        "--genetic-confounders"
            arg_type = String
            help = help = "Path to a genetic confounders file"
        "--extra-confounders"
            arg_type = String
            help = help = "Path to a confounders file to be merged with the genetic-confounders file"
        "--covariates"
            arg_type = String
            help = help = "Path to covariates file for estimation of E[Y|X]"
        "--extra-treatments"
            arg_type = String
            help = help = "Path to a additional treatments file to be merged with the genotypes file"
    end

    return parse_args(s)
end


parsed_args = parse_commandline()

finalize_tmle_inputs(parsed_args)