using ArgParse
using GenesInteraction


function parse_commandline()
    s = ArgParseSettings(
        description = "This program computes estimates of Genes Epistatis using TMLE.",
        commands_are_required = false,
        version = "0.1.0",
        add_version = true)


    @add_arg_table s begin
        "<phenotypefile>"
            help = "A PLINK phenotype file (.txt format)"
            required = true
        "<confoundersfile>"
            help = "A file (.csv format) containing the confounding variable values "*
                   " used for the analysis. The first "*
                   "line of the file should contain the columns names."
            required = true
        "<queryfile>"
            help = "A file (.csv format) containing the SNPs to consider for epistasis estimation "* 
                   "and their associated alleles that will be used as treatment and control."*
                   "The first line should contain the column names and "*
                   "each line should match the following order: \n"*
                   "SNP_1 ALLELE_TREATMENT_1 ALLELE_CONTROL_1 SNP_2 ..."
            required = true
        "<estimatorfile>"
            help = "A file (.toml format) containing the description of the tmle estimator "* 
                   "that will be used. Examples can be found in the `config` directory of the repository."
            required = true
        "<outfile>"
            help = "A path where the results will be saved (.csv format)"
            required = true
        "--verbosity", "-v"
            help = "Verbosity level"
            arg_type = Int
            default = 1
        "--imputation-threshold", "-t"
            help = "As we use imputed data, a decision must be taken regarding the most "*
                   "probable allele given the probabilities for each allele. Only most probable alleles "*
                   "passing the threshold are considered."
            arg_type = Float
            default = 0.9
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Epistasis estimation will be run with the following arguments")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    UKBBtmleepistasis(
        parsed_args["<phenotypefile>"], 
        parsed_args["<confoundersfile>"],
        parsed_args["<queryfile>"],
        parsed_args["<estimatorfile>"],
        parsed_args["<outfile>"];
        verbosity=parsed_args["verbosity"],
        threshold=parsed_args["imputation-threshold"]
    )
end

main()
