function cli_settings()
    s = ArgParseSettings(description="TargeneCore CLI.")

    @add_arg_table s begin
        "estimation-inputs"
            action = :command
            help = "Generate estimation inputs"
        
        "make-dataset"
            action = :command
            help = "Generate a dataset."

        "summary-plots"
            action = :command
            help = "Generate summary plots."
        
        "filter-chromosome"
            action = :command
            help = "Filter chromosome file based on QC filters."

        "merge-beds"
            action = :command
            help = "Merge multiple PLINK .bed files."
        
        "adapt-flashpca"
            action = :command
            help = "Adapts FlashPCA output."
    end

    @add_arg_table s["estimation-inputs"] begin
        "config-file"
            arg_type = String
            help = "Configuration file (.yaml)."
        
        "--genotypes-prefix"
            help = "Prefix path to BGEN/BED chromosomes."
            required = true
            arg_type = String

        "--traits-file"
            arg_type = String
            required = true
            help = "Path to the traits dataset"

        "--pcs-file"
            arg_type = String
            required = true
            help = "Path to a genetic confounders file"

        "--outprefix"
            arg_type = String
            help = "Output prefix"
            default = "final"

        "--batchsize"
            arg_type = Int
            help = "To split estimands in multiple files of size batchsize"

        "--call-threshold"
            arg_type = Float64
            help = "Hard call threshold for imputed genotypes"
            default = nothing

        "--positivity-constraint"
            arg_type = Float64
            default = nothing
            help = "Minimum frequency a treatment value occuring in a Parameter must reach in the population"

        "--verbosity", "-v"
            help = "Verbosity level"
            arg_type = Int
            default = 0
    end

    @add_arg_table s["make-dataset"] begin
        "--genotypes-prefix"
            help = "Prefix path to BGEN/BED chromosomes."
            required = true
            arg_type = String

        "--traits-file"
            arg_type = String
            required = true
            help = "Path to the traits dataset"

        "--pcs-file"
            arg_type = String
            required = true
            help = "Path to a genetic confounders file"

        "--variants-file"
            help = "File containing genetic variants to extract."
            required = true

        "--call-threshold"
            arg_type = Float64
            help = "Hard call threshold for imputed genotypes"
            default = 0.9

        "--out"
            help = "Output dataset file."
            arg_type = String
            default = "dataset.arrow"

        "--verbosity"
            help = "Logging level."
            arg_type = Int
            default = 0
    end

    @add_arg_table s["summary-plots"] begin
        "results-prefix"
            help = "Prefix to result files."
            required = true

        "--outprefix"
            help = "Prefix to output plots."
            default = "."

        "--verbosity"
            help = "Logging level."
            arg_type = Int
            default = 0
    end

    @add_arg_table s["filter-chromosome"] begin
        "input"
            help = "Prefix to bed files."

        "output"
            help = "Output file."
        
        "traits-file"
            arg_type = String
            required = false
            help = "Path to list of sample_ids to filter on."

        "--qc-file"
            help = "Path to the UKBiobank ukb_snp_qc.txt"
            arg_type = String
            default = nothing
            required = false

        "--maf-threshold"
            help = "SNPs with MAF lower than this value will be filtered out"
            arg_type = Float64
            default = 0.01
            required = false

        "--ld-blocks-file"
            help = "Path to the LD block definition around the SNPs of interest. "*
                   "Issued from an LD analysis."
            arg_type = String
            required = false
    end

    @add_arg_table s["merge-beds"] begin
        "input-prefix"
            help = "Prefix to bed files."

        "output"
            help = "Output file."
    end

    @add_arg_table s["adapt-flashpca"] begin
        "input"
            help = "Input file."

        "output"
            help = "Output file."
    end

    return s
end

function julia_main()::Cint
    settings = parse_args(ARGS, cli_settings())
    cmd = settings["%COMMAND%"]
    cmd_settings = settings[cmd]
    if cmd == "estimation-inputs"
        estimation_inputs(
            cmd_settings["config-file"],
            cmd_settings["genotypes-prefix"],
            cmd_settings["traits-file"],
            cmd_settings["pcs-file"],
            outprefix=cmd_settings["outprefix"],
            batchsize=cmd_settings["batchsize"],
            call_threshold=cmd_settings["call-threshold"],
            positivity_constraint=cmd_settings["positivity-constraint"],
            verbosity=cmd_settings["verbosity"]
        )
    elseif cmd == "make-dataset"
        make_dataset(
            cmd_settings["genotypes-prefix"],
            cmd_settings["traits-file"],
            cmd_settings["pcs-file"],
            cmd_settings["variants-file"];
            out=cmd_settings["out"],
            call_threshold=cmd_settings["call-threshold"],
            verbosity=cmd_settings["verbosity"]
        )
    elseif cmd == "summary-plots"
        summary_plots(
            cmd_settings["results-prefix"],
            outprefix=cmd_settings["outprefix"],
            verbosity=cmd_settings["verbosity"],
        )
    elseif cmd == "filter-chromosome"
        filter_chromosome(
            cmd_settings["input"], cmd_settings["output"], cmd_settings["traits-file"]; 
            maf_threshold=cmd_settings["maf-threshold"], 
            ld_block_file=cmd_settings["ld-blocks-file"], 
            qc_file=cmd_settings["qc-file"]
        )
    elseif cmd == "merge-beds"
        merge_beds(cmd_settings["input-prefix"], cmd_settings["output"])
    elseif cmd == "adapt-flashpca"
        adapt_flashpca(
            cmd_settings["input"], cmd_settings["output"]
        )
    else
        throw(ArgumentError(string("No function matching command:", cmd)))
    end
    return 0
end