module TestSievePlateau

using Test
using TMLEEpistasis
using DataFrames
using CSV 
using JLD2
using TMLE
using CategoricalArrays
using MLJLinearModels
using MLJBase


include("helper_fns.jl")


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

@testset "Test build_work_list" begin
    grm_ids = TMLEEpistasis.GRMIDs("data/grm/test.grm.id")
    prefix = "rs12_rs54"
    build_results_files(grm_ids, prefix)
    # mode=QUERYREPORTS, pval=0.05
    influence_curves, n_obs, file_tmlereport_pairs = TMLEEpistasis.build_work_list(prefix, grm_ids; pval=0.05)
    jldopen(string(prefix, "_batch_1_Real.hdf5")) do io
        tmlereports = io["TMLEREPORTS"]
        @test influence_curves[1, :] == convert(Vector{Float32}, tmlereports["2_1"].influence_curve)
        @test influence_curves[2, :] == convert(Vector{Float32}, tmlereports["2_2"].influence_curve)
    end
    @test n_obs == [194, 194]
    @test file_tmlereport_pairs == ["rs12_rs54_batch_1_Real.hdf5" => "2_1", "rs12_rs54_batch_1_Real.hdf5"=>"2_2"]

    # mode=MACHINE, pval=1
    influence_curves, n_obs, file_tmlereport_pairs = TMLEEpistasis.build_work_list(prefix, grm_ids; pval=1)
    @test size(influence_curves) == (6, 194)
    @test n_obs == repeat([194], 6)
    file_tmlereport_pairs == [
        "rs12_rs54_batch_1_Bool.hdf5" => "1_1",
        "rs12_rs54_batch_1_Bool.hdf5" => "1_2",
        "rs12_rs54_batch_1_Real.hdf5" => "1_1",
        "rs12_rs54_batch_1_Real.hdf5" => "1_2",
        "rs12_rs54_batch_1_Real.hdf5" => "2_1",
        "rs12_rs54_batch_1_Real.hdf5" => "2_2"
    ]
    clean_hdf5estimates_files(prefix)
end

@testset "Test bit_distance" begin
    sample_grm = Float32[-0.6, -0.8, -0.25, -0.3, -0.1, 0.1, 0.7, 0.5, 0.2, 1.]
    n??s = 6
    ??s = TMLEEpistasis.default_??s(n??s, max_??=0.75)
    @test ??s == Float32[0.0, 0.15, 0.3, 0.45, 0.6, 0.75]
    ??s = TMLEEpistasis.default_??s(n??s)
    @test ??s == Float32[0., 0.4, 0.8, 1.2, 1.6, 2.0]
    d = TMLEEpistasis.bit_distances(sample_grm, ??s)
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
    # distance indicator with 3 ??s and corresponding to row 4
    indicator = [1. 0. 0. 0.2
                 0. 0. 1. 1.
                 1. 0. 1. 1.]
    sample = 4
    var_ = TMLEEpistasis.aggregate_variances(influence_curves, indicator, sample)
    @test var_ == [24.0  189.0
                   40.0  225.0
                   48.0  333.0]
end

@testset "Test normalize!" begin
    # 2 ??s and 3 curves
    n_obs = [10, 10, 100]
    variances = [1. 2. 3.
                 4. 5. 6.]
    TMLEEpistasis.normalize!(variances, n_obs)
    @test variances == [0.1 0.2 0.03
                        0.4 0.5 0.06]
end

@testset "Test compute_variances" begin
    n_curves = 3
    n_samples = 5
    n??s = 5
    n_obs = [3, 4, 4]
    ??s = TMLEEpistasis.default_??s(n??s)
    # The GRM has 15 lower triangular elements
    grm = Float32[0.4, 0.1, 0.5, 0.2, -0.2, 0.6, 0.3, -0.6, 
                  0.4, 0.3, 0.6, 0.3, 0.7, 0.3, 0.1]
    influence_curves = Float32[0.1 0. 0.1 0.3 0.
                               0.1 0.2 0.1 0.0 0.2
                               0.0 0. 0.1 0.3 0.2]
                  
    
    variances = TMLEEpistasis.compute_variances(influence_curves, grm, ??s, n_obs)
    @test size(variances) == (n??s, n_curves)

    # when ??=2, all elements are used
    for curve_id in 1:n_curves
        s = sum(influence_curves[curve_id, :])
        var = sum(s*influence_curves[curve_id, i] for i in 1:n_samples)/n_obs[curve_id]
        @test variances[end, curve_id] ??? var
    end

    # Decreasing variances with ?? as all inf curves are positives
    for n?? in 1:n??s-1
        @test all(variances[n??, :] .<= variances[n??+1, :])
    end

    # Check against basic_variance_implementation
    matrix_distance = zeros(Float32, n_samples, n_samples)
    for ??_id in 1:n??s
        vector_distance = TMLEEpistasis.bit_distances(grm, [??s[??_id]])
        distance_vector_to_matrix!(matrix_distance, vector_distance, n_samples)
        for curve_id in 1:n_curves
            influence_curve = influence_curves[curve_id, :]
            var_ = basic_variance_implementation(matrix_distance, influence_curve, n_obs[curve_id])
            @test variances[??_id, curve_id] ??? var_
        end
    end

    # Check by hand for a single ??=0.5
    @test variances[2, :] ??? Float32[0.03666667, 0.045, 0.045]

end

@testset "Test grm_rows_bounds" begin
    n_samples = 5
    grm_bounds = TMLEEpistasis.grm_rows_bounds(n_samples)
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
    stderrors = TMLEEpistasis.corrected_stderrors(variances, n_obs)
    # sanity check
    @test size(stderrors, 1) == 10

    # check for the first curve
    stderrors[1] == sqrt(maximum(variances[:,1])/n_obs[1])

    close(io)
end

@testset "Test sieve_variance_plateau" begin
    grm_ids = TMLEEpistasis.GRMIDs("data/grm/test.grm.id")
    prefix = "rs12_rs45"
    nb_estimators = 10
    build_results_files(grm_ids, prefix)
    outfilename = "rs12_rs45_sieve_variance.hdf5"
    parsed_args = Dict(
        "prefix" => prefix,
        "pval" => 0.05,
        "grm-prefix" => "data/grm/test.grm",
        "out" => outfilename,
        "nb-estimators" => nb_estimators,
        "max-tau" => 0.75
    )

    sieve_variance_plateau(parsed_args)

    results_file = jldopen(outfilename)
    @test size(results_file["TAUS"]) == (nb_estimators,)
    @test size(results_file["VARIANCES"]) == (nb_estimators, 2)
    @test results_file["SOURCEFILE_REPORTID_PAIRS"] == 
        ["rs12_rs45_batch_1_Real.hdf5" => "2_1","rs12_rs45_batch_1_Real.hdf5" => "2_2"]
    @test size(results_file["STDERRORS"]) == (2,)

    close(results_file)

    rm(outfilename)
    clean_hdf5estimates_files(prefix)
end

end