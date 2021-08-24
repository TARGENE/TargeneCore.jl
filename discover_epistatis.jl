using ArgParse
using GenesInteraction


function parse_commandline()
    s = ArgParseSettings(
        description = "This program computes estimates of Genes Epistatis.",
        commands_are_required = false,
        version = "0.1.0",
        add_version = true)


    @add_arg_table s begin
        "<genfile>"
            help = "A PLINK Genotype file (.bed format)"
            required = true
        "<phenotypefile>"
            help = "A PLINK phenotype file (.txt format)"
            required = true
        "<confoundersfile>"
            help = "A file containing the confounding variable values "*
                   " used for the analysis (.csv format)"
            required = true
        "<snpfile>"
            help = "A file containing the SNPs to consider for interactions, "* 
                   "each line correspond to an interaction (.csv format)"
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Epistasis estimation will be run with the following arguments")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    runepistatis(
        parsed_args["<genfile>"], 
        parsed_args["<phenotypefile>"], 
        parsed_args["<confoundersfile>"],
        parsed_args["<snpfile>"];
        estimatorfile=parsed_args["<snpfile>"]
    )
end

main()
