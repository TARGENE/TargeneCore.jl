using TMLEEpistasis
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Generate queries that can be used for TMLE, see "* 
                                     "https://github.com/olivierlabayle/GenesInteraction.jl/")

    @add_arg_table s begin
        "--mode", "-m"
            help = "Specifies how the queries are generated:"*
                   "The default,`frequency` will lead to a specification of the form "*
                   "MAJOR:MAJOR -> MAJOR:MINOR for each SNP"
            arg_type = String
            default = "frequency"
        "--out", "-o"
            help = "Out directory where the queries will be saved"
            arg_type = String
            default = "."
        "--calling-threshold", "-t"
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
        "asb-prefix"
            arg_type = String
            help = "Prefix path to filtered allelic-specific binding SNPS output by baal-nf "*
                   "see https://git.ecdf.ed.ac.uk/oalmelid/baal-nf"
            required = true
        "trans-actors"
            arg_type = String
            help = "Path to trans-acting SNPS"
            required = true
        "chr-prefix"
            help = "Prefix path to BGEN chromosomes from the UKBB."
            arg_type = String
            required = true
    end

    return parse_args(s)
end


parsed_args = parse_commandline()


generate_queries(parsed_args)
