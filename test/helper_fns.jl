using BGEN

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