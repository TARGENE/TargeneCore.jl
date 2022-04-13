using TMLEEpistasis
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Generate a summary file of the full pipeline")

    @add_arg_table s begin
        "prefix"
            arg_type = String
            help = "Prefix to both estimates/sieve .hdf5 files"
            required = true
        "out"
            arg_type = String
            help = "Path where the summary file will be saved"
            required = true
        "--sieve"
            arg_type = Bool
            help = "If a sieve file is actually present in the hdf5 files"
            action = :store_true
    end

    return parse_args(s)
end


parsed_args = parse_commandline()

build_summary(parsed_args)