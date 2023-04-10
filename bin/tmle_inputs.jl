using TargeneCore
using ArgParse

function parse_commandline()
    s = ArgParseSettings(description="Preparation of data for TMLE, two modes are currently available.")

    @add_arg_table s begin
        "from-actors"
            action = :command
            help = "Will generate interaction parameters between SNPs output by the "*
                   "baal-nf pipeline (https://git.ecdf.ed.ac.uk/oalmelid/baal-nf), trans-actors "*
                   " and potential additional exposures from an external trait dataset."
        
        "from-param-files"
            action = :command
            help = "Will assume Parameter files (see README.md) to be given."

        "--out-prefix"
            arg_type = String
            help = "Output prefix"
            default = "final"

        "--traits"
            arg_type = String
            required = true
            help = "Path to the traits dataset"

        "--call-threshold"
            arg_type = Float64
            help = "This is written down on the query file, it is the threshold that"*
                    " is later used to hard call genotypes based on their probabilities"
            default = 0.9

        "--bgen-prefix"
            help = "Prefix path to BGEN chromosomes."
            required = true
            arg_type = String

        "--pcs"
            arg_type = String
            required = true
            help = "Path to a genetic confounders file"

        "--phenotype-batch-size"
            arg_type = Int
            required = false
            help = "Further performance may be obtained by batching phenotypes in a"*
                   " single Targeted Estimation run. If not specified, all phenotypes "*
                   "constitute a batch."

        "--positivity-constraint"
            arg_type = Float64
            default = 0.01
            required = false
            help = "Minimum frequency a treatment value occuring in a Parameter must reach in the population"
    end

    @add_arg_table s["from-actors"] begin
        "bqtls"
            arg_type = String
            help = "Path to bqtls"

        "trans-actors-prefix"
            arg_type = String
            help = "Prefix path to trans-acting SNPs files"
        
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

        "--genotypes-as-int"
            action = :store_true
            help = "If true, genotypes are encoded as the number of minor alleles."
    end

    @add_arg_table s["from-param-files"] begin
        "param-prefix"
            arg_type = String
            help = "Prefix to parameter files"
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

tmle_inputs(parsed_args)