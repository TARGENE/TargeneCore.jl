module TestTMLEInputs

using Test
using CSV
using DataFrames
using TargeneCore
using YAML
using StableRNGs
using BGEN

function cleanup()
    for file in readdir()
        if startswith(file, "final.")
            rm(file)
        end
    end
end

#####################################################################
##################           UNIT TESTS            ##################
#####################################################################

@testset "Test asb_snps / trans_actors" begin
    bQTLs = TargeneCore.asb_snps(joinpath("data", "asb_files", "asb_"))
    @test bQTLs == ["RSID_17", "RSID_99", "RSID_198"]

    eQTLs = TargeneCore.trans_actors(joinpath("data", "trans_actors_fake.csv"))
    @test eQTLs == ["RSID_102", "RSID_2"]
end

@testset "Test genotypes_encoding" begin
    b = Bgen(BGEN.datadir("example.8bits.bgen"))
    v = variant_by_rsid(b, "RSID_10")
    minor_allele_dosage!(b, v)
    # The minor allele is the first one
    @test minor_allele(v) == alleles(v)[1]
    @test TargeneCore.genotypes_encoding(v) == [2, 1, 0]

    # The minor allele is the second one
    v = variant_by_rsid(b, "RSID_102")
    minor_allele_dosage!(b, v)
    @test minor_allele(v) == alleles(v)[2]
    @test TargeneCore.genotypes_encoding(v) == [0, 1, 2]
end

@testset "Test call_genotypes for a single SNP" begin
    probabilities = [NaN 0.3 0.2 0.9;
                     NaN 0.5 0.2 0.05;
                     NaN 0.2 0.6 0.05]
    variant_genotypes = [2, 1, 0]

    threshold = 0.9
    genotypes = TargeneCore.call_genotypes(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1] === genotypes[2] === genotypes[3] === missing
    @test genotypes[4] == 2

    threshold = 0.55
    genotypes = TargeneCore.call_genotypes(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1] === genotypes[2]  === missing
    @test genotypes[3] == 0
    @test genotypes[4] == 2
end

@testset "Test call_genotypes for all SNPs" begin
    snp_list = ["RSID_10", "RSID_100"]
    genotypes = TargeneCore.call_genotypes(joinpath("data", "ukbb", "imputed" , "ukbb"), snp_list, 0.95)
    # I only look at the first 10 rows
    # SAMPLE_ID    
    @test genotypes[1:9, "SAMPLE_ID"] == ["sample_00$i" for i in 1:9]
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == ones(10)
    # RSID_100
    @test all(genotypes[1:10, "RSID_100"] .=== [1, 2, 1, missing, 1, 1, missing, 1, 0, 1])
    # Test column order
    @test DataFrames.names(genotypes) == ["SAMPLE_ID", "RSID_10", "RSID_100"]
end

@testset "Test build_asb_trans_param_files!" begin
    rng = StableRNG(123)
    param_files = Pair{String, Dict}[]
    n = 100
    treatments = DataFrame(
        RSID_1 = rand(rng, [2, 1, 0], n),
        RSID_2 = rand(rng, [0, 1, 2], n),
        EXTRA_TREAT = rand(rng, [0,1], n)
    )
    tuple_1 = ["RSID_1", "RSID_2", "EXTRA_TREAT"]
    freqs = TargeneCore.frequency_table(treatments, tuple_1)
    @test freqs == Dict{Any, Any}(
        (0, 0, 0) => 0.09,
        (1, 2, 1) => 0.07,
        (0, 2, 1) => 0.1,
        (1, 1, 0) => 0.07,
        (0, 1, 0) => 0.04,
        (1, 2, 0) => 0.05,
        (2, 0, 1) => 0.02,
        (0, 2, 0) => 0.05,
        (2, 1, 1) => 0.05,
        (2, 0, 0) => 0.03,
        (2, 2, 1) => 0.06,
        (2, 1, 0) => 0.04,
        (1, 0, 1) => 0.06,
        (0, 0, 1) => 0.05,
        (1, 1, 1) => 0.06,
        (0, 1, 1) => 0.05,
        (2, 2, 0) => 0.06,
        (1, 0, 0) => 0.05)
    treatment_tuple_generator = [tuple_1]
    TargeneCore.build_asb_trans_param_files!(param_files, treatment_tuple_generator, treatments; positivity_constraint=0.)
    filename, param_file = only(param_files)
    @test filename == "I_RSID_1_RSID_2_EXTRA_TREAT.yaml"
    @test param_file["Treatments"] == tuple_1
    # There are 3 * 3 * 1 expected parameters
    @test param_file["Parameters"] == [
        Dict("name" => "I__RSID_1_0->1__RSID_2_0->1__EXTRA_TREAT_0->1",
        "RSID_1" => Dict("case" => 1, "control" => 0), 
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 1, "control" => 0)),
        Dict("name" => "I__RSID_1_0->2__RSID_2_0->1__EXTRA_TREAT_0->1",
        "RSID_1" => Dict("case" => 2, "control" => 0), 
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 1, "control" => 0)),
        Dict("name" => "I__RSID_1_1->2__RSID_2_0->1__EXTRA_TREAT_0->1", 
        "RSID_1" => Dict("case" => 2, "control" => 1), 
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 1, "control" => 0)),
        Dict("name" => "I__RSID_1_0->1__RSID_2_0->2__EXTRA_TREAT_0->1", 
        "RSID_1" => Dict("case" => 1, "control" => 0), 
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 2, "control" => 0)),
        Dict("name" => "I__RSID_1_0->2__RSID_2_0->2__EXTRA_TREAT_0->1",
        "RSID_1" => Dict("case" => 2, "control" => 0),  
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 2, "control" => 0)),
        Dict("name" => "I__RSID_1_1->2__RSID_2_0->2__EXTRA_TREAT_0->1", 
        "RSID_1" => Dict("case" => 2, "control" => 1), 
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 2, "control" => 0)),
        Dict("name" => "I__RSID_1_0->1__RSID_2_1->2__EXTRA_TREAT_0->1", 
        "RSID_1" => Dict("case" => 1, "control" => 0), 
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 2, "control" => 1)),
        Dict("name" => "I__RSID_1_0->2__RSID_2_1->2__EXTRA_TREAT_0->1", 
        "RSID_1" => Dict("case" => 2, "control" => 0), 
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 2, "control" => 1)),
        Dict("name" => "I__RSID_1_1->2__RSID_2_1->2__EXTRA_TREAT_0->1", 
        "RSID_1" => Dict("case" => 2, "control" => 1), 
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 2, "control" => 1))
    ]
    param_files = Pair{String, Dict}[]
    TargeneCore.build_asb_trans_param_files!(param_files, treatment_tuple_generator, treatments; positivity_constraint=0.04)
    filename, param_file = only(param_files)
    @test filename == "I_RSID_1_RSID_2_EXTRA_TREAT.yaml"
    @test param_file["Treatments"] == tuple_1
    @test param_file["Parameters"] == [
        Dict("name" => "I__RSID_1_0->1__RSID_2_0->1__EXTRA_TREAT_0->1", 
        "RSID_1" => Dict("case" => 1, "control" => 0),     
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 1, "control" => 0)),
        Dict("name" => "I__RSID_1_0->1__RSID_2_0->2__EXTRA_TREAT_0->1", 
        "RSID_1" => Dict("case" => 1, "control" => 0), 
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 2, "control" => 0)),
        Dict("name" => "I__RSID_1_0->1__RSID_2_1->2__EXTRA_TREAT_0->1", 
        "RSID_1" => Dict("case" => 1, "control" => 0), 
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 2, "control" => 1)),
        Dict("name" => "I__RSID_1_0->2__RSID_2_1->2__EXTRA_TREAT_0->1", 
        "RSID_1" => Dict("case" => 2, "control" => 0), 
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0), 
        "RSID_2" => Dict("case" => 2, "control" => 1)),
        Dict("name" => "I__RSID_1_1->2__RSID_2_1->2__EXTRA_TREAT_0->1",
        "RSID_1" => Dict("case" => 2, "control" => 1),  
        "EXTRA_TREAT" => Dict("case" => 1, "control" => 0),
        "RSID_2" => Dict("case" => 2, "control" => 1))
    ]
end

@testset "Test treatments_from_actors" begin
    # At least two type of actors should be specified
    @test_throws ArgumentError TargeneCore.treatments_from_actors(nothing, nothing, 1)
    @test_throws ArgumentError TargeneCore.treatments_from_actors(1, nothing, nothing)
    @test_throws ArgumentError TargeneCore.treatments_from_actors(nothing, 1, nothing)

    bqtl_file = joinpath("data", "bqtls.csv")
    trans_actors_prefix = joinpath("data", "trans_actors_1.csv")
    env_file = joinpath("data", "env.txt")
    # bqtls and trans_actors
    bqtls, transactors, envs = TargeneCore.treatments_from_actors(bqtl_file, nothing, trans_actors_prefix)
    @test bqtls isa DataFrame
    @test envs isa Nothing
    @test transactors isa Vector{DataFrame}
    @test size(transactors, 1) == 1

    # bqtls and env
    bqtls, transactors, envs = TargeneCore.treatments_from_actors(bqtl_file, env_file, nothing)
    @test bqtls isa DataFrame
    @test envs == DataFrame(ENV = ["sex"])
    @test transactors isa Nothing

    # trans actors and env
    bqtls, transactors, envs = TargeneCore.treatments_from_actors(nothing, env_file, trans_actors_prefix)
    @test bqtls isa Nothing
    @test envs == DataFrame(ENV = ["sex"])
    @test transactors isa Vector{DataFrame}
    @test size(transactors, 1) == 1
end

@testset "Test combine_by_bqtl" begin
    bQTLS = DataFrame(
        ID = ["rs1", "rs2", "rs3"],
        CHR = [1, 2, 3],
        DESCRIPTION = ["VDR", "VDR", "VDR"]
    )
    trans_actors = [
        DataFrame(
            ID = ["rs4", "rs5"],
            CHR = [1, 12],
            DESCRIPTION = ["VitD-QTL", "VitD-QTL"]
    ),
        DataFrame(
            ID = ["rs6", "rs7", "rs8"],
            CHR = [12, 9, 8],
            DESCRIPTION = ["RXR-QTL", "RXR-QTL", "RXR-QTL"]
        )
    ]
    envs = DataFrame(ID=["sex"], DESCRIPTION=["Env"])
    # order 2, no environmental variable
    combinations = TargeneCore.combine_by_bqtl(bQTLS, trans_actors, nothing, 2)
    @test combinations == [Dict("T" => ["rs1", "rs4"]),
                           Dict("T" => ["rs1", "rs5"]),
                           Dict("T" => ["rs2", "rs4"]),
                           Dict("T" => ["rs2", "rs5"]),
                           Dict("T" => ["rs3", "rs4"]),
                           Dict("T" => ["rs3", "rs5"]),
                           Dict("T" => ["rs1", "rs6"]),
                           Dict("T" => ["rs1", "rs7"]),
                           Dict("T" => ["rs1", "rs8"]),
                           Dict("T" => ["rs2", "rs6"]),
                           Dict("T" => ["rs2", "rs7"]),
                           Dict("T" => ["rs2", "rs8"]),
                           Dict("T" => ["rs3", "rs6"]),
                           Dict("T" => ["rs3", "rs7"]),
                           Dict("T" => ["rs3", "rs8"])]
    # order 2, only one trans-actors dataset
    combinations = TargeneCore.combine_by_bqtl(bQTLS, [trans_actors[1]], nothing, 2)
    @test combinations == [Dict("T" => ["rs1", "rs4"]),
                           Dict("T" => ["rs1", "rs5"]),
                           Dict("T" => ["rs2", "rs4"]),
                           Dict("T" => ["rs2", "rs5"]),
                           Dict("T" => ["rs3", "rs4"]),
                           Dict("T" => ["rs3", "rs5"])]
    # order 3, no environmental variable
    combinations = TargeneCore.combine_by_bqtl(bQTLS, trans_actors, nothing, 3)
    @test combinations == [Dict("T" => ["rs1", "rs4", "rs6"]),
                           Dict("T" => ["rs1", "rs4", "rs7"]),
                           Dict("T" => ["rs1", "rs4", "rs8"]),
                           Dict("T" => ["rs1", "rs5", "rs6"]),
                           Dict("T" => ["rs1", "rs5", "rs7"]),
                           Dict("T" => ["rs1", "rs5", "rs8"]),
                           Dict("T" => ["rs2", "rs4", "rs6"]),
                           Dict("T" => ["rs2", "rs4", "rs7"]),
                           Dict("T" => ["rs2", "rs4", "rs8"]),
                           Dict("T" => ["rs2", "rs5", "rs6"]),
                           Dict("T" => ["rs2", "rs5", "rs7"]),
                           Dict("T" => ["rs2", "rs5", "rs8"]),
                           Dict("T" => ["rs3", "rs4", "rs6"]),
                           Dict("T" => ["rs3", "rs4", "rs7"]),
                           Dict("T" => ["rs3", "rs4", "rs8"]),
                           Dict("T" => ["rs3", "rs5", "rs6"]),
                           Dict("T" => ["rs3", "rs5", "rs7"]),
                           Dict("T" => ["rs3", "rs5", "rs8"])]
    
    # order 2, environmental variable
    combinations = TargeneCore.combine_by_bqtl(bQTLS, trans_actors, envs, 2)
    @test combinations == [Dict("T" => ["rs1", "rs4"]),
                           Dict("T" => ["rs1", "rs5"]),
                           Dict("T" => ["rs2", "rs4"]),
                           Dict("T" => ["rs2", "rs5"]),
                           Dict("T" => ["rs3", "rs4"]),
                           Dict("T" => ["rs3", "rs5"]),
                           Dict("T" => ["rs1", "rs6"]),
                           Dict("T" => ["rs1", "rs7"]),
                           Dict("T" => ["rs1", "rs8"]),
                           Dict("T" => ["rs2", "rs6"]),
                           Dict("T" => ["rs2", "rs7"]),
                           Dict("T" => ["rs2", "rs8"]),
                           Dict("T" => ["rs3", "rs6"]),
                           Dict("T" => ["rs3", "rs7"]),
                           Dict("T" => ["rs3", "rs8"]),
                           Dict("T" => ["rs1", "sex"]),
                           Dict("T" => ["rs2", "sex"]),
                           Dict("T" => ["rs3", "sex"])]
    # order 3, only one trans-actors dataset
    combinations = TargeneCore.combine_by_bqtl(bQTLS, [trans_actors[1]], envs, 3)
    @test combinations == [Dict("T" => ["rs1", "rs4", "sex"]),
                           Dict("T" => ["rs1", "rs5", "sex"]),
                           Dict("T" => ["rs2", "rs4", "sex"]),
                           Dict("T" => ["rs2", "rs5", "sex"]),
                           Dict("T" => ["rs3", "rs4", "sex"]),
                           Dict("T" => ["rs3", "rs5", "sex"])]
    # order 2, no trans-actors dataset
    combinations = TargeneCore.combine_by_bqtl(bQTLS, nothing, envs, 2)
    @test combinations == [Dict("T" => ["rs1", "sex"]),
                           Dict("T" => ["rs2", "sex"]),
                           Dict("T" => ["rs3", "sex"])]
end
#####################################################################
###############           END-TO-END TESTS            ###############
#####################################################################

@testset "Test tmle_inputs with-param-files: scenario 1" begin
    # Scenario:
    # - binary and continuous phenotypes
    # - genetic and extra confounders
    # - no covariates
    # - no extra treatments
    # - no batch size
    # - no positivity constraint
    parsed_args = Dict(
        "with-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_")), 
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => nothing, 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => joinpath("data", "extra_confounders.csv"), 
        "%COMMAND%" => "with-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "genetic-confounders" => joinpath("data", "genetic_confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => nothing,
        "phenotype-batch-size" => nothing,
        "positivity-constraint" => 0.
    )
    @test_logs(
        (:warn, "Some treatment variables could not be read from the data files and associated parameter files will not be processed: TREAT_1"), 
        tmle_inputs(parsed_args)
    )

    # Data Files
    confounders = CSV.read("final.confounders.csv", DataFrame)
    @test names(confounders) == ["SAMPLE_ID", "PC1", "PC2", "21003", "22001"]
    @test size(confounders) == (490, 5)

    treatments = CSV.read("final.treatments.csv", DataFrame)
    @test size(treatments) == (490, 3)
    @test names(treatments) == ["SAMPLE_ID", "RSID_2", "RSID_198"]
    @test Set(unique(treatments[!, "RSID_2"])) == Set([missing, 0, 1, 2])
    @test Set(unique(treatments[!, "RSID_198"])) == Set([0, 1, 2])

    binary_phenotypes = CSV.read("final.binary-phenotypes.csv", DataFrame)
    @test size(binary_phenotypes) == (490, 3)
    @test names(binary_phenotypes) == ["SAMPLE_ID", "BINARY_1", "BINARY_2"]

    continuous_phenotypes = CSV.read("final.continuous-phenotypes.csv", DataFrame)  
    @test size(continuous_phenotypes) == (490, 3)
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "CONTINUOUS_1", "CONTINUOUS_2"]

    # Parameter files: untouched because no phenotype batches were specified
    # They are deduplicated for both continuous and binary phenotypes
    @test YAML.load_file(joinpath("config", "param_1.yaml")) == YAML.load_file("final.binary.param_1.yaml") == YAML.load_file("final.continuous.param_1.yaml")
    @test YAML.load_file(joinpath("config", "param_2.yaml")) == YAML.load_file("final.binary.param_2.yaml") == YAML.load_file("final.continuous.param_2.yaml")
    
    cleanup()
end

@testset "Test tmle_inputs with-param-files: scenario 2" begin
    # Scenario:
    # - binary and continuous phenotypes
    # - no extra confounders
    # - covariates
    # - SNP and extra treatments
    # - batch size: 1
    # - no positivity constraint
    parsed_args = Dict(
        "with-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_1")), 
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => joinpath("data", "extra_treatments.csv"), 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => nothing, 
        "%COMMAND%" => "with-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "genetic-confounders" => joinpath("data", "genetic_confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => joinpath("data", "covariates.csv"),
        "phenotype-batch-size" => 1,
        "positivity-constraint" => 0.
    )
    tmle_inputs(parsed_args)

    confounders = CSV.read("final.confounders.csv", DataFrame)
    @test names(confounders) == ["SAMPLE_ID", "PC1", "PC2"]
    @test size(confounders) == (490, 3)

    treatments = CSV.read("final.treatments.csv", DataFrame)
    @test size(treatments) == (490, 3)
    @test names(treatments) == ["SAMPLE_ID", "RSID_2", "TREAT_1"]

    binary_phenotypes = CSV.read("final.binary-phenotypes.csv", DataFrame)
    @test size(binary_phenotypes) == (490, 3)
    @test names(binary_phenotypes) == ["SAMPLE_ID", "BINARY_1", "BINARY_2"]

    continuous_phenotypes = CSV.read("final.continuous-phenotypes.csv", DataFrame)  
    @test size(continuous_phenotypes) == (490, 3)
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "CONTINUOUS_1", "CONTINUOUS_2"]

    covariates = CSV.read("final.covariates.csv", DataFrame)
    @test names(covariates) == ["SAMPLE_ID", "COV_1"]
    @test size(covariates) == (490, 2)
    # Parameter files: modified because phenotype batches was specified
    # They are deduplicated for both continuous and binary phenotypes
    origin_1 = YAML.load_file(joinpath("config", "param_1.yaml"))
    binary_1_1 = YAML.load_file("final.binary.param_1.batch_1.yaml")
    binary_1_2 = YAML.load_file("final.binary.param_1.batch_2.yaml")
    continuous_1_1 = YAML.load_file("final.continuous.param_1.batch_1.yaml")
    continuous_1_2 = YAML.load_file("final.continuous.param_1.batch_2.yaml")
    # Parameters and Treatments sections unchanged
    @test binary_1_1["Parameters"] == binary_1_2["Parameters"] == continuous_1_1["Parameters"] == origin_1["Parameters"]
    @test binary_1_1["Treatments"] == binary_1_2["Treatments"] == continuous_1_2["Treatments"] == origin_1["Treatments"]
    # Targets sections changed
    @test binary_1_1["Targets"] == ["BINARY_1"]
    @test binary_1_2["Targets"] == ["BINARY_2"]
    @test continuous_1_1["Targets"] == ["CONTINUOUS_1"]
    @test continuous_1_2["Targets"] == ["CONTINUOUS_2"]

    origin_2 = YAML.load_file(joinpath("config", "param_1_with_extra_treatment.yaml"))
    binary_2_1 = YAML.load_file("final.binary.param_1_with_extra_treatment.batch_1.yaml")
    binary_2_2 = YAML.load_file("final.binary.param_1_with_extra_treatment.batch_2.yaml")
    continuous_2_1 = YAML.load_file("final.continuous.param_1_with_extra_treatment.batch_1.yaml")
    continuous_2_2 = YAML.load_file("final.continuous.param_1_with_extra_treatment.batch_2.yaml")
    # Parameters and Treatments sections unchanged
    @test binary_2_1["Parameters"] == binary_2_2["Parameters"] == continuous_2_1["Parameters"] == origin_2["Parameters"]
    @test binary_2_1["Treatments"] == binary_2_2["Treatments"] == continuous_2_2["Treatments"] == origin_2["Treatments"]
    # Targets sections changed
    @test binary_2_1["Targets"] == ["BINARY_1"]
    @test binary_2_2["Targets"] == ["BINARY_2"]
    @test continuous_2_1["Targets"] == ["CONTINUOUS_1"]
    @test continuous_2_2["Targets"] == ["CONTINUOUS_2"]

    cleanup()
end

@testset "Test tmle_inputs with-asb-trans: scenario 1" begin
    # Scenario:
    # - binary and continuous phenotypes
    # - genetic and extra confounders
    # - no covariates
    # - no extra treatments
    # - no batch size
    # - no param
    # - no positivity constraint
    parsed_args = Dict(
        "with-asb-trans" => Dict{String, Any}(
            "asb-prefix" => joinpath("data", "asb_files", "asb"), 
            "trans-actors" => joinpath("data", "trans_actors_fake.csv"),
            "param-prefix" => nothing
            ),
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => nothing, 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => joinpath("data", "extra_confounders.csv"), 
        "%COMMAND%" => "with-asb-trans", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "genetic-confounders" => joinpath("data", "genetic_confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => nothing,
        "phenotype-batch-size" => nothing,
        "positivity-constraint" => 0.
    )
    tmle_inputs(parsed_args)

    confounders = CSV.read("final.confounders.csv", DataFrame)
    @test names(confounders) == ["SAMPLE_ID", "PC1", "PC2", "21003", "22001"]
    @test size(confounders) == (490, 5)

    treatments = CSV.read("final.treatments.csv", DataFrame)
    @test size(treatments) == (490, 6)
    @test names(treatments) == ["SAMPLE_ID", "RSID_2", "RSID_102", "RSID_17", "RSID_198", "RSID_99"]

    binary_phenotypes = CSV.read("final.binary-phenotypes.csv", DataFrame)
    @test size(binary_phenotypes) == (490, 3)
    @test names(binary_phenotypes) == ["SAMPLE_ID", "BINARY_1", "BINARY_2"]

    continuous_phenotypes = CSV.read("final.continuous-phenotypes.csv", DataFrame)  
    @test size(continuous_phenotypes) == (490, 3)
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "CONTINUOUS_1", "CONTINUOUS_2"]

    snp_combinations = Set(Iterators.product(["RSID_102", "RSID_2"], ["RSID_17", "RSID_198", "RSID_99"]))
    for snp_combination in copy(snp_combinations)
        file_pattern = string("I_", join(snp_combination, "_"), ".yaml")
        binary_file = YAML.load_file(string("final.binary.", file_pattern))
        continuous_file = YAML.load_file(string("final.continuous.", file_pattern))
        @test binary_file == continuous_file
        eqtl, bqtl = Tuple(binary_file["Treatments"])
        setdiff!(snp_combinations, [(eqtl, bqtl)])
        if bqtl == "RSID_198"
            @test size(binary_file["Parameters"], 1) == 5
        else
            @test size(binary_file["Parameters"], 1) == 3
        end
    end
    @test snp_combinations == Set()
    cleanup()
end

@testset "Test tmle_inputs with-asb-trans: scenario 2" begin
    # Scenario:
    # - binary and continuous phenotypes
    # - no extra confounders
    # - covariates
    # - extra treatments
    # - batch size
    # - param
    # - no positivity constraint
    parsed_args = Dict(
        "with-asb-trans" => Dict{String, Any}(
            "asb-prefix" => joinpath("data", "asb_files", "asb"), 
            "trans-actors" => joinpath("data", "trans_actors_fake.csv"),
            "param-prefix" => joinpath("config", "template")
            ),
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => joinpath("data", "extra_treatments.csv"), 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => nothing, 
        "%COMMAND%" => "with-asb-trans", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "genetic-confounders" => joinpath("data", "genetic_confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => joinpath("data", "covariates.csv"),
        "phenotype-batch-size" => 1,
        "positivity-constraint" => 0.0
    )
    tmle_inputs(parsed_args)

    confounders = CSV.read("final.confounders.csv", DataFrame)
    @test names(confounders) == ["SAMPLE_ID", "PC1", "PC2"]
    @test size(confounders) == (490, 3)
    
    treatments = CSV.read("final.treatments.csv", DataFrame)
    @test size(treatments) == (490, 7)
    @test names(treatments) == ["SAMPLE_ID", "RSID_2", "RSID_102", "RSID_17", "RSID_198", "RSID_99", "TREAT_1"]
    
    binary_phenotypes = CSV.read("final.binary-phenotypes.csv", DataFrame)
    @test size(binary_phenotypes) == (490, 3)
    @test names(binary_phenotypes) == ["SAMPLE_ID", "BINARY_1", "BINARY_2"]
    
    continuous_phenotypes = CSV.read("final.continuous-phenotypes.csv", DataFrame)  
    @test size(continuous_phenotypes) == (490, 3)
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "CONTINUOUS_1", "CONTINUOUS_2"]
    
    covariates = CSV.read("final.covariates.csv", DataFrame)
    @test names(covariates) == ["SAMPLE_ID", "COV_1"]
    @test size(covariates) == (490, 2)
    
    # There are two parameter files, one with extra treatments and one without
    # For each phenotype type (e.g. continuous or binary):
    # We thus expect a maximum of 2 * nb_phenotypes * n_eqtls * n_bqtls = 24 parameter files 
    treatment_combinations = Set(Iterators.product(["RSID_102", "RSID_2"], ["RSID_17", "RSID_198", "RSID_99"], ["TREAT_1"]))
    for treatment_combination in copy(treatment_combinations)
        file_pattern = string("I_", join(treatment_combination, "_"))
        binary_file_1 = YAML.load_file(string("final.binary.", file_pattern, ".batch_1.yaml"))
        binary_file_2 = YAML.load_file(string("final.binary.", file_pattern, ".batch_2.yaml"))
        continuous_file_1 = YAML.load_file(string("final.continuous.", file_pattern, ".batch_1.yaml"))
        continuous_file_2 = YAML.load_file(string("final.continuous.", file_pattern, ".batch_2.yaml"))
        @test binary_file_1["Parameters"] == continuous_file_1["Parameters"]
        @test binary_file_1["Treatments"] == continuous_file_1["Treatments"]
        @test binary_file_2["Parameters"] == continuous_file_2["Parameters"]
        @test binary_file_2["Treatments"] == continuous_file_2["Treatments"]
        setdiff!(treatment_combinations, [treatment_combination])
    end
    @test treatment_combinations == Set()
    cleanup()

    # Adding positivity constraint, only 20 files are generated
    parsed_args["positivity-constraint"] = 0.01
    tmle_inputs(parsed_args)

    param_files = filter(x -> occursin.(r"^final.*yaml$",x), readdir())
    @test size(param_files, 1) == 40 

    cleanup()
end




end

true