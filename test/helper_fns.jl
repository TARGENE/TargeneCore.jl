using BGEN
using Random
using MLJLinearModels
using CategoricalArrays
using TMLE
using JLD2

###############################################################################
# JLD2Saver Callback (Dup from TargetedEstimation otherwise can't merge projects)
###############################################################################

"""
A callback to save the TMLE estimation results to disk instead of keeping them in memory
"""
mutable struct JLD2Saver <: TMLE.Callback
    file::String
    save_machines::Bool
    save_ic::Bool
end

"""
    maybe_summarize(report::TMLEReport, save_ic)

This is a temporary way to remove the influence curve data to save disk space
while waiting for the new TMLE API.
"""
function maybe_summarize(report::TMLEReport, save_ic)
    if save_ic
        return report
    else
        return TMLE.summarize(report)
    end
end

function TMLE.after_tmle(callback::JLD2Saver, report::TMLEReport, target_id::Int, query_id::Int)
    jldopen(callback.file, "a"; compress=true) do io
        group = haskey(io, "TMLEREPORTS") ? io["TMLEREPORTS"] : JLD2.Group(io, "TMLEREPORTS")
        group[string(target_id, "_", query_id)] = maybe_summarize(report, callback.save_ic)
    end
end

function TMLE.after_fit(callback::JLD2Saver, mach::TMLE.Machine, id::Symbol)
    if callback.save_machines
        jldopen(callback.file, "a"; compress=true) do io
            group = haskey(io, "MACHINES") ? io["MACHINES"] : JLD2.Group(io, "MACHINES")
            group[string(id)] = mach
        end
    end
end

function TMLE.finalize(callback::JLD2Saver, estimation_report::NamedTuple)
    jldopen(callback.file, "a"; compress=true) do io
        io["low_propensity_scores"] = estimation_report.low_propensity_scores
    end
    return estimation_report
end

function clean_hdf5estimates_files(prefix)
    rm(string(prefix, "_batch_1_Real.hdf5"))
    rm(string(prefix, "_batch_1_Bool.hdf5"))
    rm(string(prefix, "_batch_2_Bool.hdf5"))
end


function build_results_files(grm_ids, prefix; save_ic=true)
    rng = Xoshiro(0)
    n = size(grm_ids, 1)
    T = (t₁=categorical(rand(rng, [0, 1], n)),)
    W = (w₁=rand(rng, n), w₂=rand(rng, n))
    height = rand(rng, n)
    bmi = 2convert(Vector{Float64}, T.t₁) + 0.2rand(rng, n)
    cancer = categorical(rand(rng, [0,1], n))
    query_1 = Query(case=(t₁=1,), control=(t₁=0,), name="QUERY_1")
    query_2 = Query(case=(t₁=0,), control=(t₁=1,), name="QUERY_2")

    # Continuous single batch
    cont_path = string(prefix, "_batch_1_Real.hdf5")
    tmle_reg = TMLEstimator(
            LinearRegressor(), 
            LogisticClassifier(),
            query_1,
            query_2
    )
    TMLE.fit(tmle_reg, T, W, (height=height, bmi=bmi);
        verbosity=0,
        cache=false,
        callbacks=[JLD2Saver(cont_path, true, save_ic)])
    if save_ic
        jldopen(cont_path, "a") do io
                io["SAMPLE_IDS"] = Dict(
                    "height" => string.(grm_ids.SAMPLE_ID),
                    "bmi" => string.(grm_ids.SAMPLE_ID)
                )
        end
    end
    # Binary single batch
    bin_path = string(prefix, "_batch_1_Bool.hdf5")
    tmle_bin = TMLEstimator(
        LogisticClassifier(), 
        LogisticClassifier(),
        query_1,
        query_2)
    TMLE.fit(tmle_bin, T, W, (cancer=cancer,);
        verbosity=0,
        cache=false,
        callbacks=[JLD2Saver(bin_path, true, save_ic)])
    if save_ic
        jldopen(bin_path, "a") do io
            io["SAMPLE_IDS"] = Dict(
                "cancer" => string.(grm_ids.SAMPLE_ID)
            )
        end
    end
    # Empty file that may be output by the tmle procedure
    jldopen(string(prefix, "_batch_2_Bool.hdf5"), "a") do io
        JLD2.Group(io, "TMLEREPORTS")
    end
end


function test_queries(queries, expected_queries)
    for (i, query) in enumerate(queries)
        @test query.case == expected_queries[i].case
        @test query.control == expected_queries[i].control
        @test query.name == expected_queries[i].name
    end
end


function _test_wiped_data(mach)
    @test mach.data == ()
    @test :data ∉ mach.cache
    @test mach.args == ()
end


function test_data_has_been_removed(tmle_mach)
    _test_wiped_data(tmle_mach)

    for submach in tmle_mach.report.machines
        _test_wiped_data(submach)
    end

    for submach in report(tmle_mach).G.machines
        _test_wiped_data(submach)
    end

    for submach in report(tmle_mach).Q̅.machines
        _test_wiped_data(submach)
    end
end