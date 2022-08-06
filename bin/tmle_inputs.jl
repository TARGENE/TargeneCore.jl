using TargeneCore
using ArgParse

function parse_commandline()
    s = ArgParseSettings(description="Preparation of data for TMLE, two modes are currently available.")

    @add_arg_table s begin
        "with-asb-trans"
            action = :command
            help = "Will generate interaction parameters between SNPs output by the "*
                   "baal-nf pipeline (https://git.ecdf.ed.ac.uk/oalmelid/baal-nf), trans-actors "*
                   " and potential additional exposures from an external trait dataset."
        
        "with-param-files"
            action = :command
            help = "Will assume Parameter files (see README.md) to be given."

        "--out-prefix"
            arg_type = String
            help = "Output prefix"
            default = "final"

        "--binary-phenotypes"
            arg_type = String
            help = "Path to the binary phenotypes file"

        "--continuous-phenotypes"
            arg_type = String
            help = "Path to the continuous phenotypes file"
        
        "--call-threshold"
            arg_type = Float64
            help = "This is written down on the query file, it is the threshold that"*
                    " is later used to hard call genotypes based on their probabilities"
            default = 0.9

        "--bgen-prefix"
            help = "Prefix path to BGEN chromosomes."
            arg_type = String

        "--genetic-confounders"
            arg_type = String
            help = "Path to a genetic confounders file"

        "--extra-confounders"
            arg_type = String
            help = "Path to a confounders file to be merged with the genetic-confounders file"

        "--covariates"
            arg_type = String
            help = "Path to covariates file for estimation of E[Y|X]"

        "--extra-treatments"
            arg_type = String
            help = "Path to a additional treatments file to be merged with the genotypes file"
        
        "--phenotype-batch-size"
            arg_type = Int
            required = false
            help = "Further performance may be obtained by batching phenotypes in a"*
                   " single Targeted Estimation run"
        "--positivity-constraint"
            arg_type = Float64
            default = 0.01
            required = false
            help = "Minimum frequency a treatment value occuring in a Parameter must reach in the population"

    end

    @add_arg_table s["with-asb-trans"] begin
        "asb-prefix"
            arg_type = String
            help = "Prefix path to filtered allelic-specific binding SNPS output by baal-nf "*
                   "see https://git.ecdf.ed.ac.uk/oalmelid/baal-nf"
        "trans-actors"
            arg_type = String
            help = "Path to trans-acting SNPS"
        "--param-prefix"
            arg_type = String
            help = "Prefix to template parameter file to provide additional treatment variable information"
            required = false
    end

    @add_arg_table s["with-param-files"] begin
        "--param-prefix"
            arg_type = String
            help = "Prefix to parameter files"
            required = true
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

finalize_tmle_inputs(parsed_args)