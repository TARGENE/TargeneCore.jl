using ArgParse
using GenesInteraction


function parse_commandline()
    s = ArgParseSettings(
        description = "This program computes estimates of Genes Epistatis.",
        commands_are_required = false,
        version = "0.1.0",
        add_version = true)


    @add_arg_table s begin
        "<snpfile>"
            help = "The SNP file .bed"
            required = true
        "<phenotypefile>"
            help = "The phenotype file"
            required = true
        "<confoundersfile>"
            help = "The confounders file"
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

    runepistatis(parsed_args["<snpfile>"], parsed_args["<phenotypefile>"], parsed_args["<confoundersfile>"])
end

main()
