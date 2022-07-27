using Test
using TMLEEpistasis

function read_output(parsed_args)
    treatments = CSV.read(
        string(parsed_args["out-prefix"], ".treatments.csv"),
        DataFrame
    )
    phenotypes = CSV.read(
        string(parsed_args["out-prefix"], ".phenotypes.csv"),
        DataFrame
    )
    confounders = CSV.read(
        string(parsed_args["out-prefix"], ".confounders.csv"),
        DataFrame
    )
    return treatments, phenotypes, confounders
end

@testset "Test finalize_tmle_inputs" begin
    parsed_args = Dict{String, Union{String, Nothing}}(
        "out-prefix" => "final",
        "phenotypes" => "test.phenotypes.csv",
        "genotypes" => "test.genotypes.csv",
        "covariates" => "test.covariates.csv",
        "extra-confounders" => "test.extra-confounders.csv",
        "genetic-confounders" => "test.genetic-confounders.csv",
        "extra-treatments" => "test.extra-treatments.csv"
    )
    # Write Data
    CSV.write(parsed_args["phenotypes"], 
        DataFrame(
            SAMPLE_ID = [1, 2, 3, 4, 5],
            PHENO_1 = [1, 2, 3, 4, 5],
            PHENO_2 = [1, 1, 1, 1, 1]
    ))
    CSV.write(parsed_args["genetic-confounders"], 
        DataFrame(
            SAMPLE_ID = [2, 3, 4, 5, 1],
            PC_1 = [0.2, 0.4, 0.1, 0.3, 0.8]
    ))
    CSV.write(parsed_args["extra-confounders"], 
        DataFrame(
            SAMPLE_ID = [5, 3, 4, 2, 1],
            CONF_1 = [0.2, 0.4, 0.1, 0.3, 0.8]
    ))
    CSV.write(parsed_args["covariates"], 
        DataFrame(
            SAMPLE_ID = [1, 3, 5, 2, 4],
            COV_1 = [0.2, 0.4, 0.1, 0.3, 0.8]
    ))
    CSV.write(parsed_args["extra-treatments"], 
        DataFrame(
            SAMPLE_ID = [3, 5, 2, 1, 4],
            T_1 = [0.2, 0.4, 0.1, 0.3, 0.8]
    ))
    CSV.write(parsed_args["genotypes"], 
        DataFrame(
            SAMPLE_ID = [4, 5, 2, 1, 3],
            G_1 = ["AC", "CC", "CC", "AA", "AA"]
    ))
    # First scenario
    parsed_args["extra-confounders"] = nothing
    parsed_args["extra-treatments"] = nothing
    parsed_args["covariates"] = nothing
    finalize_tmle_inputs(parsed_args)
    treatments, phenotypes, confounders = read_output(parsed_args)

    @test treatments.SAMPLE_ID == phenotypes.SAMPLE_ID == confounders.SAMPLE_ID
    @test names(phenotypes) == ["SAMPLE_ID", "PHENO_1", "PHENO_2"]
    @test names(treatments) == ["SAMPLE_ID", "G_1"]
    @test names(confounders) == ["SAMPLE_ID", "PC_1"]
    @test !isfile("final.covariates.csv")

    # Second scenario
    parsed_args["extra-confounders"] = "test.extra-confounders.csv"
    parsed_args["extra-treatments"] = "test.extra-treatments.csv"
    finalize_tmle_inputs(parsed_args)
    treatments, phenotypes, confounders = read_output(parsed_args)

    @test treatments.SAMPLE_ID == phenotypes.SAMPLE_ID == confounders.SAMPLE_ID
    @test names(phenotypes) == ["SAMPLE_ID", "PHENO_1", "PHENO_2"]
    @test names(treatments) == ["SAMPLE_ID", "G_1", "T_1"]
    @test names(confounders) == ["SAMPLE_ID", "PC_1", "CONF_1"]
    @test !isfile("final.covariates.csv")

    # Third scenario
    parsed_args["covariates"] = "test.covariates.csv"
    finalize_tmle_inputs(parsed_args)
    treatments, phenotypes, confounders = read_output(parsed_args)
    covariates  = CSV.read(
        string(parsed_args["out-prefix"], ".covariates.csv"),
        DataFrame
    )
    
    @test treatments.SAMPLE_ID == phenotypes.SAMPLE_ID == confounders.SAMPLE_ID == covariates.SAMPLE_ID
    @test names(phenotypes) == ["SAMPLE_ID", "PHENO_1", "PHENO_2"]
    @test names(treatments) == ["SAMPLE_ID", "G_1", "T_1"]
    @test names(confounders) == ["SAMPLE_ID", "PC_1", "CONF_1"]
    @test names(covariates) == ["SAMPLE_ID", "COV_1"]

    # Remove input files
    for (key, filepath) in parsed_args
        if key != "out-prefix"
            rm(filepath)
        end
    end
    # Remove output files
    for key in ("confounders", "covariates", "phenotypes", "treatments")
        rm(string(parsed_args["out-prefix"], ".", key, ".csv"))
    end
end