using ArgParse
using TargeneCore

function parse_commandline()
    s = ArgParseSettings(
        description = "Summary plot generation.",
        commands_are_required = false,
        version = "0.1",
        add_version = true)

    @add_arg_table s begin
        "results-prefix"
            help = "Prefix to result files."
            required = true

        "--out-prefix"
            help = "Prefix to output plots."
            default = "."

        "--verbosity"
            help = "Logging level."
            arg_type = Int
            default = 1
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

generate_summary_plots(parsed_args)

