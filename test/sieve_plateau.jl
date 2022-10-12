module TestSievePlateau

using Test
using TargeneCore
using DataFrames
using CSV 
using JLD2
using TMLE
using CategoricalArrays
using TargetedEstimation
using StableRNGs
using Distributions
using LogExpFunctions

function build_dataset(sample_ids)
    rng = StableRNG(123)
    n = size(sample_ids, 1)
    # Confounders
    W₁ = rand(rng, Uniform(), n)
    W₂ = rand(rng, Uniform(), n)
    # Covariates
    C₁ = rand(rng, n)
    # Treatment | Confounders
    T₁ = rand(rng, Uniform(), n) .< logistic.(0.5sin.(W₁) .- 1.5W₂)
    T₂ = rand(rng, Uniform(), n) .< logistic.(-3W₁ - 1.5W₂)
    # target | Confounders, Covariates, Treatments
    μ = 1 .+ 2W₁ .+ 3W₂ .- 4C₁.*T₁ .+ T₁ + T₂.*W₂.*T₁
    y₁ = μ .+ rand(rng, Normal(0, 0.01), n)
    y₂ = rand(rng, Uniform(), n) .< logistic.(μ)
    # Add some missingness
    y₂ = vcat(missing, y₂[2:end])

    dataset = DataFrame(
        SAMPLE_ID = string.(sample_ids),
        T1 = categorical(T₁),
        T2 = categorical(T₂),
        W1 = W₁, 
        W2 = W₂,
        C1 = C₁,
        CONTINUOUS_TARGET = y₁,
        BINARY_TARGET = categorical(y₂)
    )

    CSV.write("data.csv", dataset)
end

function build_tmle_output_file(sample_ids, param_file, outfile)
    build_dataset(sample_ids)
    # Only one continuous phenotype / machines not saved / no adaptive cv
    parsed_args = Dict(
        "data" => "data.csv",
        "param-file" => param_file,
        "estimator-file" => joinpath("config", "tmle_config.yaml"),
        "out" => outfile,
        "verbosity" => 0,
        "save-full" => true,
    )

    TargetedEstimation.main(parsed_args)
end


function basic_variance_implementation(matrix_distance, influence_curve, n_obs)
    variance = 0.f0
    n_samples = size(influence_curve, 1)
    for i in 1:n_samples
        for j in 1:n_samples
            variance += matrix_distance[i, j]*influence_curve[i]* influence_curve[j]
        end
    end
    variance/n_obs
end

function distance_vector_to_matrix!(matrix_distance, vector_distance, n_samples)
    index = 1
    for i in 1:n_samples
        for j in 1:i
            # enforce indicator = 1 when i =j 
            if i == j
                matrix_distance[i, j] = 1
            else
                matrix_distance[i, j] = vector_distance[index]
                matrix_distance[j, i] = vector_distance[index]
            end
            index += 1
        end
    end
end

function test_initial_output(output, expected_output)
    # Metadata columns
    for col in [:PARAMETER_TYPE, :TREATMENTS, :CASE, :CONTROL, :TARGET, :CONFOUNDERS, :COVARIATES]
        for index in eachindex(output[!, col])
            if expected_output[index, col] === missing
                @test expected_output[index, col] === output[index, col]
            else
                @test expected_output[index, col] == output[index, col]
            end
        end
    end
    # TMLE i.i.d columns
    for col in (:INITIAL_ESTIMATE, :ESTIMATE, :STD, :PVALUE, :LWB, :UPB)
        @test eltype(output[!, col]) == Float64
    end
    # Sieve columns
    for col in (:SIEVE_STD, :SIEVE_PVALUE, :SIEVE_LWB, :SIEVE_UPB)
        @test all(x===missing for x in output[!, col])
    end
end


@testset "Test readGRM" begin
    prefix = joinpath("data", "grm", "test.grm")
    GRM, ids = TargeneCore.readGRM(prefix)
    @test eltype(ids.SAMPLE_ID) == String
    @test size(GRM, 1) == 18915
    @test size(ids, 1) == 194
end

@testset "Test build_work_list" begin
    grm_ids = TargeneCore.GRMIDs(joinpath("data", "grm", "test.grm.id"))
    param_file_1 = joinpath("config", "sieve_tests_parameters_1.yaml")
    outfile_1 = "tmle_output_1.hdf5"
    prefix = "tmle_output"
    # CASE_1: only one file
    build_tmle_output_file(grm_ids.SAMPLE_ID, param_file_1, outfile_1)
    # Since pval = 1., all parameters are considered for sieve variance
    output, influence_curves, n_obs, row_ids = TargeneCore.build_work_list(prefix, grm_ids; pval=1.)
    @test row_ids == [1, 2, 3, 4, 5, 6]
    @test n_obs == [194, 194, 194, 193, 193, 193]
    # Check influence curves
    io = jldopen(outfile_1)
    continuous_results = io["results"]["CONTINUOUS_TARGET"]["tmle_results"]
    for i in 1:3
        @test convert(Vector{Float32}, continuous_results[i].IC) == influence_curves[i, :]
    end
    binary_results = io["results"]["BINARY_TARGET"]["tmle_results"]
    for i in 1:3
        @test convert(Vector{Float32}, vcat(0, binary_results[i].IC)) == influence_curves[i+3, :]
    end
    # Check output
    some_expected_cols = DataFrame(
        PARAMETER_TYPE = ["IATE", "IATE", "ATE", "IATE", "IATE", "ATE"],
        TREATMENTS = ["T2_&_T1", "T2_&_T1", "T2_&_T1", "T2_&_T1", "T2_&_T1", "T2_&_T1"],
        CASE = ["1_&_1", "0_&_1", "1_&_1", "1_&_1", "0_&_1", "1_&_1"],
        CONTROL = ["0_&_0", "1_&_0", "0_&_0", "0_&_0", "1_&_0", "0_&_0"],
        TARGET = ["CONTINUOUS_TARGET", "CONTINUOUS_TARGET", "CONTINUOUS_TARGET", "BINARY_TARGET", "BINARY_TARGET", "BINARY_TARGET"],
        CONFOUNDERS = ["W1_&_W2", "W1_&_W2", "W1_&_W2", "W1_&_W2", "W1_&_W2", "W1_&_W2"],
        COVARIATES = ["C1", "C1", "C1", "C1", "C1", "C1"]
    )
    test_initial_output(output, some_expected_cols)

    # CASE_2: add another file 
    param_file_2 = joinpath("config", "sieve_tests_parameters_2.yaml")
    outfile_2 = "tmle_output_2.hdf5"
    build_tmle_output_file(grm_ids.SAMPLE_ID, param_file_2, outfile_2)
    # This p-value filters the influence curves for the binary outcome
    output, influence_curves, n_obs, row_ids = TargeneCore.build_work_list(prefix, grm_ids; pval=1e-10)
    @test size(influence_curves) == (5, 194)
    @test n_obs == [194, 194, 194, 194, 194]
    @test row_ids == [1, 2, 3, 7, 8]

    # Check output
    some_expected_cols = DataFrame(
        PARAMETER_TYPE = ["IATE", "IATE", "ATE", "IATE", "IATE", "ATE", "ATE", "CM"],
        TREATMENTS = ["T2_&_T1", "T2_&_T1", "T2_&_T1", "T2_&_T1", "T2_&_T1", "T2_&_T1", "T1", "T1"],
        CASE = ["1_&_1", "0_&_1", "1_&_1", "1_&_1", "0_&_1", "1_&_1", "1", "0"],
        CONTROL = ["0_&_0", "1_&_0", "0_&_0", "0_&_0", "1_&_0", "0_&_0", "0", missing],
        TARGET = ["CONTINUOUS_TARGET", "CONTINUOUS_TARGET", "CONTINUOUS_TARGET", "BINARY_TARGET", "BINARY_TARGET", "BINARY_TARGET", "CONTINUOUS_TARGET", "CONTINUOUS_TARGET"],
        CONFOUNDERS = ["W1_&_W2", "W1_&_W2", "W1_&_W2", "W1_&_W2", "W1_&_W2", "W1_&_W2", "W1", "W1"],
        COVARIATES = ["C1", "C1", "C1", "C1", "C1", "C1", missing, missing]
    )
    test_initial_output(output, some_expected_cols)

    # CASE_3: empty influence curves
    output, influence_curves, n_obs, row_ids = TargeneCore.build_work_list(prefix, grm_ids; pval=0.)
    @test influence_curves == []
    @test row_ids == []
    @test n_obs == []
    test_initial_output(output, some_expected_cols)
    # clean
    rm(outfile_1)
    rm(outfile_2)
    rm("data.csv")
end

@testset "Test bit_distance" begin
    sample_grm = Float32[-0.6, -0.8, -0.25, -0.3, -0.1, 0.1, 0.7, 0.5, 0.2, 1.]
    nτs = 6
    τs = TargeneCore.default_τs(nτs, max_τ=0.75)
    @test τs == Float32[0.0, 0.15, 0.3, 0.45, 0.6, 0.75]
    τs = TargeneCore.default_τs(nτs)
    @test τs == Float32[0., 0.4, 0.8, 1.2, 1.6, 2.0]
    d = TargeneCore.bit_distances(sample_grm, τs)
    @test d == [0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0
                0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0
                0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0
                0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0
                1.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
                1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0]
end

@testset "Test aggregate_variances" begin
    # 2 influence curves containing 5 individuals
    influence_curves = [1. 2. 3. 4. 5.
                        6. 7. 8. 9. 10.]
    # distance indicator with 3 τs and corresponding to row 4
    indicator = [1. 0. 0. 0.2
                 0. 0. 1. 1.
                 1. 0. 1. 1.]
    sample = 4
    var_ = TargeneCore.aggregate_variances(influence_curves, indicator, sample)
    @test var_ == [24.0  189.0
                   40.0  225.0
                   48.0  333.0]
end

@testset "Test normalize!" begin
    # 2 τs and 3 curves
    n_obs = [10, 10, 100]
    variances = [1. 2. 3.
                 4. 5. 6.]
    TargeneCore.normalize!(variances, n_obs)
    @test variances == [0.1 0.2 0.03
                        0.4 0.5 0.06]
end

@testset "Test compute_variances" begin
    n_curves = 3
    n_samples = 5
    nτs = 5
    n_obs = [3, 4, 4]
    τs = TargeneCore.default_τs(nτs)
    # The GRM has 15 lower triangular elements
    grm = Float32[0.4, 0.1, 0.5, 0.2, -0.2, 0.6, 0.3, -0.6, 
                  0.4, 0.3, 0.6, 0.3, 0.7, 0.3, 0.1]
    influence_curves = Float32[0.1 0. 0.1 0.3 0.
                               0.1 0.2 0.1 0.0 0.2
                               0.0 0. 0.1 0.3 0.2]
                  
    
    variances = TargeneCore.compute_variances(influence_curves, grm, τs, n_obs)
    @test size(variances) == (nτs, n_curves)

    # when τ=2, all elements are used
    for curve_id in 1:n_curves
        s = sum(influence_curves[curve_id, :])
        var = sum(s*influence_curves[curve_id, i] for i in 1:n_samples)/n_obs[curve_id]
        @test variances[end, curve_id] ≈ var
    end

    # Decreasing variances with τ as all inf curves are positives
    for nτ in 1:nτs-1
        @test all(variances[nτ, :] .<= variances[nτ+1, :])
    end

    # Check against basic_variance_implementation
    matrix_distance = zeros(Float32, n_samples, n_samples)
    for τ_id in 1:nτs
        vector_distance = TargeneCore.bit_distances(grm, [τs[τ_id]])
        distance_vector_to_matrix!(matrix_distance, vector_distance, n_samples)
        for curve_id in 1:n_curves
            influence_curve = influence_curves[curve_id, :]
            var_ = basic_variance_implementation(matrix_distance, influence_curve, n_obs[curve_id])
            @test variances[τ_id, curve_id] ≈ var_
        end
    end

    # Check by hand for a single τ=0.5
    @test variances[2, :] ≈ Float32[0.03666667, 0.045, 0.045]

end

@testset "Test grm_rows_bounds" begin
    n_samples = 5
    grm_bounds = TargeneCore.grm_rows_bounds(n_samples)
    @test grm_bounds == [1 => 1
                         2 => 3
                         4 => 6
                         7 => 10
                         11 => 15]
end

@testset "Test corrected_stderrors" begin
    io = jldopen(joinpath("data", "sieve_variances.hdf5"))
    variances = io["variances"]
    n_obs = [10, 10, 10, 10, 10, 100, 100, 1000, 1000, 1000]
    stderrors = TargeneCore.corrected_stderrors(variances, n_obs)
    # sanity check
    @test size(stderrors, 1) == 10

    # check for the first curve
    stderrors[1] == sqrt(maximum(variances[:,1])/n_obs[1])

    close(io)
end

@testset "Test sieve_variance_plateau" begin
    # Generate data
    nb_estimators = 10
    grm_ids = TargeneCore.GRMIDs(joinpath("data", "grm", "test.grm.id"))
    param_file_1 = joinpath("config", "sieve_tests_parameters_1.yaml")
    tmle_outfile_1 = "tmle_output_1.hdf5"
    param_file_2 = joinpath("config", "sieve_tests_parameters_2.yaml")
    tmle_outfile_2 = "tmle_output_2.hdf5"
    build_tmle_output_file(grm_ids.SAMPLE_ID, param_file_1, tmle_outfile_1)
    build_tmle_output_file(grm_ids.SAMPLE_ID, param_file_2, tmle_outfile_2)

    outprefix = "sieve_output"
    parsed_args = Dict(
        "prefix" => "tmle_output",
        "pval" => 1e-10,
        "grm-prefix" => "data/grm/test.grm",
        "out-prefix" => outprefix,
        "nb-estimators" => nb_estimators,
        "max-tau" => 0.75,
        "verbosity" => 0
    )

    sieve_variance_plateau(parsed_args)
    # check hdf5 file
    io = jldopen(string(outprefix, ".hdf5"))
    @test io["taus"] == TargeneCore.default_τs(nb_estimators; max_τ=parsed_args["max-tau"])
    @test size(io["variances"]) == (10, 5)
    close(io)
    # check csv file
    output = CSV.read(string(outprefix, ".csv"), DataFrame)
    for col in (:SIEVE_PVALUE, :SIEVE_STD, :SIEVE_LWB, :SIEVE_UPB)
        @test all(x===missing for x in output[[4, 5, 6], col])
        @test all(x isa Float64 for x in output[[1, 2, 3, 7, 8], col])
    end

    # clean
    rm(string(outprefix, ".csv"))
    rm(string(outprefix, ".hdf5"))
    rm(tmle_outfile_1)
    rm(tmle_outfile_2)
    rm("data.csv")
end

@testset "Test sieve_variance_plateau: empty work list" begin
    # Generate data
    nb_estimators = 10
    grm_ids = TargeneCore.GRMIDs(joinpath("data", "grm", "test.grm.id"))
    param_file_1 = joinpath("config", "sieve_tests_parameters_1.yaml")
    tmle_outfile_1 = "tmle_output_1.hdf5"
    build_tmle_output_file(grm_ids.SAMPLE_ID, param_file_1, tmle_outfile_1)

    outprefix = "sieve_output"
    parsed_args = Dict(
        "prefix" => "tmle_output",
        "pval" => 0.,
        "grm-prefix" => "data/grm/test.grm",
        "out-prefix" => outprefix,
        "nb-estimators" => nb_estimators,
        "max-tau" => 0.75,
        "verbosity" => 0
    )

    sieve_variance_plateau(parsed_args)
    # check hdf5 file
    io = jldopen(string(outprefix, ".hdf5"))
    @test io["taus"] == TargeneCore.default_τs(nb_estimators; max_τ=parsed_args["max-tau"])
    close(io)
    # check csv file
    output = CSV.read(string(outprefix, ".csv"), DataFrame)
    for col in (:SIEVE_PVALUE, :SIEVE_STD, :SIEVE_LWB, :SIEVE_UPB)
        @test all(x===missing for x in output[!, col])
    end

    # clean
    rm(string(outprefix, ".csv"))
    rm(string(outprefix, ".hdf5"))
    rm(tmle_outfile_1)
    rm("data.csv")
end

end