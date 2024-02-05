using TargeneCore
using ArgParse

function parse_commandline()
    s = ArgParseSettings(description="Preparation of data for TMLE, two modes are currently available.")

    @add_arg_table s begin
        "from-actors"
            action = :command
            help = "Will generate interaction estimands between SNPs output by the "*
                   "baal-nf pipeline (https://git.ecdf.ed.ac.uk/oalmelid/baal-nf), trans-actors "*
                   " and potential additional exposures from an external trait dataset."
        
        "from-param-file"
            action = :command
            help = "Will assume Parameter files (see README.md) to be given."

        "allele-independent"
            action = :command
            help = "Generates allele independent estimands from a configuration file."
        
        "permutation-tests"
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

    @add_arg_table s["from-actors"] begin
        "bqtls"
            arg_type = String
            help = "Path to bqtls"

        "trans-actors-prefix"
            arg_type = String
            help = "Prefix path to trans-acting SNPs files"
            "--traits"
            arg_type = String
            required = true
            help = "Path to the traits dataset"

        "--call-threshold"
            arg_type = Float64
            help = "Hard call threshold for imputed genotypes"
            default = 0.9

        "--bgen-prefix"
            help = "Prefix path to BGEN chromosomes."
            required = true
            arg_type = String

        "--pcs"
            arg_type = String
            required = true
            help = "Path to a genetic confounders file"

        "--extra-confounders"
            arg_type = String
            required = false
            help = "Path to a confounders file"

        "--extra-covariates"
            arg_type = String
            required = false
            help = "Path to extra covariates file"

        "--extra-treatments"
            arg_type = String
            required = false
            help = "Path to a additional treatments file"
        
        "--orders"
            arg_type = String
            required = false
            help = "Interaction orders to be estimated"
            default = "1,2"
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
            default = 0.9

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
            default = 0.9

        "--bgen-prefix"
            help = "Prefix path to BGEN chromosomes."
            required = true
            arg_type = String

        "--pcs"
            arg_type = String
            required = true
            help = "Path to a genetic confounders file"
    end

    @add_arg_table s["permutation-tests"] begin
        "dataset"
            help = "Path to tmle input dataset file (.csv|.arrow)"
            required = true
        "results"
            help = "CSV file containing TarGene results."
            required = true
        "--pval-threshold"
            help = "The p-value threshold for significant results calling"
            default = 0.05
            arg_type = Float64
        "--estimator-key"
            help = "Estimator to use to check significance."
            default = "TMLE"
            arg_type = String
        "--limit"
            help = "The max number of permutation parameters to be generated"
            default = nothing
            arg_type = Int
        "--rng"
            help = "The random seed for permutations"
            default = 123
            arg_type = Int
        "--orders"
            help = "A comma separated set of combination orders e.g. 1,2,3"
            default = "1"
            arg_type = String
        "--max-attempts"
            help = "Maximum number of permutation to induce to maximize the number of permuted dataset matching the positivity constraint"
            default = 1
            arg_type = Int
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

tl_inputs(parsed_args)