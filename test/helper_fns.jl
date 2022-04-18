using BGEN
using Random
using MLJBase
using MLJLinearModels
using TMLE
using JLD2

genotypesfile = joinpath("data", "genotypes.csv")
confoundersfile = joinpath("data", "confounders.csv")
queryfile = joinpath("config", "query.toml")
continuous_phenotypefile = joinpath("data", "continuous_phenotypes.csv")
binary_phenotypefile = joinpath("data", "binary_phenotypes.csv")
phenotypelist_file = joinpath("data", "phen_list.csv")


function build_query_file(;path=queryfile, threshold=0.9, het="AG")
    queries = Dict(
        "threshold" => threshold,
        "SNPS" => Dict(
            "RSID_10"  => BGEN.datadir("example.8bits.bgen"),
            "RSID_100" => BGEN.datadir("example.8bits.bgen")),
        "QUERY_1" => Dict(
            "RSID_10"  => "GG -> "*het,
            "RSID_100" => "GG -> "*het),
        "QUERY_2" => Dict(
            "RSID_10"  => "GG -> "*het,
            "RSID_100" => "GG -> AA"
        )
    )
    open(path, "w") do io
        TOML.print(io, queries)
    end
end


function build_ate_query_file(;path=queryfile, threshold=0.9)
    queries = Dict(
        "threshold" => threshold,
        "SNPS" => Dict("RSID_10"  => BGEN.datadir("example.8bits.bgen")),
        "QUERY_1" => Dict("RSID_10"  => "GG -> AG"),
        "QUERY_2" => Dict("RSID_10"  => "GG -> AA")
    )
    open(path, "w") do io
        TOML.print(io, queries)
    end
end


function clean_hdf5estimates_files(prefix)
    rm(string(prefix, "_batch_1_Real.hdf5"))
    rm(string(prefix, "_batch_1_Bool.hdf5"))
end


function build_results_files(grm_ids, prefix; mode="QUERYREPORTS")
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
    tmle_reg = TMLEstimator(
        LinearRegressor(), 
        LogisticClassifier(),
        query_1,
        query_2)
    mach_reg = machine(tmle_reg, T, W, (height=height, bmi=bmi), cache=false)
    fit!(mach_reg, verbosity=0)
    cont_path = string(prefix, "_batch_1_Real.hdf5")
    jldopen(cont_path, "w") do io
        io["SAMPLE_IDS"] = Dict(
            "height" => string.(grm_ids.SAMPLE_ID),
            "bmi" => string.(grm_ids.SAMPLE_ID)
        )
        if mode == "MACHINE"
            io["MACHINE"] = mach_reg
        else
            io["QUERYREPORTS"] = queryreports(mach_reg)
        end
    end

    # Binary single batch
    tmle_bin = TMLEstimator(
        LogisticClassifier(), 
        LogisticClassifier(),
        query_1,
        query_2)
    mach_bin = machine(tmle_bin, T, W, (cancer=cancer,), cache=false)
    fit!(mach_bin, verbosity=0)
    bin_path = string(prefix, "_batch_1_Bool.hdf5")
    jldopen(bin_path, "w") do io
        io["SAMPLE_IDS"] = Dict(
            "cancer" => string.(grm_ids.SAMPLE_ID)
        )
        if mode == "MACHINE"
            io["MACHINE"] = mach_bin
        else
            io["QUERYREPORTS"] = queryreports(mach_bin)
        end
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