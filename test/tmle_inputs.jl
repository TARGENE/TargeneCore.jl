using Test
using CSV
using DataFrames
using TargeneCore

function read_output(parsed_args)
    treatments = CSV.read(
        string(parsed_args["out-prefix"], ".treatments.csv"),
        DataFrame
    )
    binary_phenotypes = CSV.read(
        string(parsed_args["out-prefix"], ".binary-phenotypes.csv"),
        DataFrame
    )
    continuous_phenotypes = CSV.read(
        string(parsed_args["out-prefix"], ".continuous-phenotypes.csv"),
        DataFrame
    )
    confounders = CSV.read(
        string(parsed_args["out-prefix"], ".confounders.csv"),
        DataFrame
    )
    return treatments, binary_phenotypes, continuous_phenotypes, confounders
end

@testset "Test tmle_inputs with-param-files" begin
    parsed_args = Dict(
        "exclude" => nothing, 
        "with-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_")), 
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => nothing, 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => nothing, 
        "%COMMAND%" => "with-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "minor-genotype-freq" => 0.001, 
        "genetic-confounders" => joinpath("data", "confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => nothing,
        "phenotype-batch-size" => nothing
    )
    tmle_inputs(parsed_args)
end

@testset "TO change" begin

    parsed_args = Dict{Any, Any}(
        "out-prefix" => "final",
        "binary-phenotypes" => "test.binary-phenotypes.csv",
        "continuous-phenotypes" => "test.continuous-phenotypes.csv",
        "genotypes" => "test.genotypes.csv",
        "covariates" => "test.covariates.csv",
        "extra-confounders" => "test.extra-confounders.csv",
        "genetic-confounders" => "test.genetic-confounders.csv",
        "extra-treatments" => "test.extra-treatments.csv"
    )
    # Write Data
    CSV.write(parsed_args["binary-phenotypes"], 
        DataFrame(
            SAMPLE_ID = [1, 2, 3, 4, 5],
            PHENO_1 = [1, 2, 3, 4, 5],
    ))
    CSV.write(parsed_args["continuous-phenotypes"], 
    DataFrame(
        SAMPLE_ID = [1, 2, 3, 4, 5],
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
    treatments, binary_phenotypes, continuous_phenotypes, confounders = read_output(parsed_args)

    @test treatments.SAMPLE_ID == binary_phenotypes.SAMPLE_ID == continuous_phenotypes.SAMPLE_ID == confounders.SAMPLE_ID
    @test names(binary_phenotypes) == ["SAMPLE_ID", "PHENO_1"]
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "PHENO_2"]
    @test names(treatments) == ["SAMPLE_ID", "G_1"]
    @test names(confounders) == ["SAMPLE_ID", "PC_1"]
    @test !isfile("final.covariates.csv")

    # Second scenario
    parsed_args["extra-confounders"] = "test.extra-confounders.csv"
    parsed_args["extra-treatments"] = "test.extra-treatments.csv"
    finalize_tmle_inputs(parsed_args)
    treatments, binary_phenotypes, continuous_phenotypes, confounders = read_output(parsed_args)

    @test treatments.SAMPLE_ID == binary_phenotypes.SAMPLE_ID == continuous_phenotypes.SAMPLE_ID == confounders.SAMPLE_ID
    @test names(binary_phenotypes) == ["SAMPLE_ID", "PHENO_1"]
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "PHENO_2"]
    @test names(treatments) == ["SAMPLE_ID", "G_1", "T_1"]
    @test names(confounders) == ["SAMPLE_ID", "PC_1", "CONF_1"]
    @test !isfile("final.covariates.csv")

    # Third scenario
    parsed_args["covariates"] = "test.covariates.csv"
    finalize_tmle_inputs(parsed_args)
    treatments, binary_phenotypes, continuous_phenotypes, confounders = read_output(parsed_args)
    covariates  = CSV.read(
        string(parsed_args["out-prefix"], ".covariates.csv"),
        DataFrame
    )
    
    @test treatments.SAMPLE_ID == binary_phenotypes.SAMPLE_ID == continuous_phenotypes.SAMPLE_ID == confounders.SAMPLE_ID == covariates.SAMPLE_ID
    @test names(binary_phenotypes) == ["SAMPLE_ID", "PHENO_1"]
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "PHENO_2"]
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
    for key in ("confounders", "covariates", "binary-phenotypes", "continuous-phenotypes","treatments")
        rm(string(parsed_args["out-prefix"], ".", key, ".csv"))
    end
end