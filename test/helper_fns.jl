
confoundersfile = joinpath("data", "confounders.csv")
queryfile = joinpath("config", "query.toml")
phenotypefile = joinpath("data", "phenotypes.csv")
phenotypelist_file_1 = joinpath("data", "phen_list_1.csv")
phenotypelist_file_2 = joinpath("data", "phen_list_2.csv")

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
