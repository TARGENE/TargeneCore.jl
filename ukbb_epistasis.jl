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
        "phenotypes"
            help = "A file (.csv format). The first row contains the column names with `eid` the sample ID"
            required = true
        "confounders"
            help = "A file (.csv format) containing the confounding variables values and the sample ids associated"*
                   " with the participants. The first line of the file should contain the columns names and the sample ids "*
                   " column name should be: `SAMPLE_ID`."
            required = true
        "queries"
            help = "A file (.toml format) see: config/sample_query.toml for more information"
            required = true
        "estimator"
            help = "A file (.toml format) describing the tmle estimator to use, see config/test_categorical.toml"*
                   " or config/test_continuous.toml for a basic example."
            required = true
        "output"
            help = "A path where the results will be saved (.csv format)"
            required = true
        "--phenotype", "-p"
            help = "The phenotype to consider for the analysis"
            required = true
            arg_type = String
        "--verbosity", "-v"
            help = "Verbosity level"
            arg_type = Int
            default = 1
    end

    return parse_args(s)
end

parsed_args = parse_commandline()


TMLEEpistasisUKBB(parsed_args)