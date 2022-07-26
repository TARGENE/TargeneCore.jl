using TMLEEpistasis
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Generate queries that can be used for TMLE, see "* 
                                     "https://github.com/olivierlabayle/GenesInteraction.jl/")

    @add_arg_table s begin
        "--mode", "-m"
            help = "`given` or `asb`. If `given` then the option `--query-prefix` must be"*
                   " provided. If `asb` options `--asb-prefix` and `--trans-actors must be`"*
                   " provided. Queries will be generated as MAJOR:MAJOR -> MAJOR:MINOR for each SNP."
            arg_type = String
            default = "given"
        "--outdir", "-o"
            help = "Out directory where the queries will be saved"
            arg_type = String
            default = "."
        "--query-prefix", "-q"
            help = "Prefix path to queries, must be provided if mode is: `given`"
            arg_type = String
            required = false
        "--call-threshold", "-c"
            arg_type = Float64
            help = "This is written down on the query file, it is the threshold that"*
                    " is later used to hard call genotypes based on their probabilities"
            default = 0.9
        "--exclude", "-e"
            arg_type = String
            help = "A file containing SNPs to be excluded from the analysis (one per line)"
            required = false
        "--minor-genotype-freq", "-f"
            arg_type = Float64
            help = "Minor genotype frequency filter"
            default = 0.001
            required = false
        "--asb-prefix", "-a"
            arg_type = String
            help = "Prefix path to filtered allelic-specific binding SNPS output by baal-nf "*
                   "see https://git.ecdf.ed.ac.uk/oalmelid/baal-nf"
            required = false
        "--trans-actors", "-t"
            arg_type = String
            help = "Path to trans-acting SNPS"
            required = false
        "--sample_ids"
            arg_type = String
            help = "Path to list of sample_ids to filter on."
            required = false
        "chr-prefix"
            help = "Prefix path to BGEN chromosomes from the UKBB."
            arg_type = String
            required = true
    end

    return parse_args(s)
end


parsed_args = parse_commandline()


build_genotypes_and_queries(parsed_args)
