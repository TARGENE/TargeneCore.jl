using ArgParse
using GenesInteraction


function parse_commandline()
    s = ArgParseSettings(
        description = "This program computes estimates of Genes Epistatis from the UK-Biobank using TMLE."*
                      "Here is a list of the arguments that should be provided, you can also have a look at the "*
                      "test/data and test/config folders to see some examples.",
        commands_are_required = false,
        version = "0.1.0",
        add_version = true)


    @add_arg_table s begin
        "<phenotypefile>"
            help = "A file (.csv format) as output by the `ukbconv` program. The first trait appearing in the file will be used."*
                   "It is assumed that the sample ids are given by a column named `eid`"
            required = true
        "<confoundersfile>"
            help = "A file (.csv format) containing the confounding variables values and the sample ids associated"*
                   " with the participants. The first line of the file should contain the columns names and the sample ids "*
                   " column name should be: `SAMPLE_ID`."
            required = true
        "<queryfile>"
            help = "A file (.toml format) see: config/sample_query.toml for more information"
            required = true
        "<estimatorfile>"
            help = "A file (.toml format) describing the tmle estimator to use, see config/test_categorical.toml"*
                   " or config/test_continuous.toml for a basic example."
            required = true
        "<outfile>"
            help = "A path where the results will be saved (.csv format)"
            required = true
        "--verbosity", "-v"
            help = "Verbosity level"
            arg_type = Int
            default = 1
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Epistasis estimation will be run with the following arguments")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    TMLEEpistasisUKBB(
        parsed_args["<phenotypefile>"], 
        parsed_args["<confoundersfile>"],
        parsed_args["<queryfile>"],
        parsed_args["<estimatorfile>"],
        parsed_args["<outfile>"];
        verbosity=parsed_args["verbosity"],
    )
end

main()
