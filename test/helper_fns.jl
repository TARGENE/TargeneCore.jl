
confoundersfile = joinpath("data", "confounders.csv")
queryfile = joinpath("config", "query.toml")
phenotypefile = joinpath("data", "phenotypes.csv")

function build_query_file(;path=queryfile, threshold=0.9, het="AG")
    queries = Dict(
        "threshold" => threshold,
        "SNPS" => Dict(
            "RSID_10"  => BGEN.datadir("example.8bits.bgen"),
            "RSID_100" => BGEN.datadir("example.8bits.bgen")),
        "QUERY_1" => Dict(
            "RSID_10"  => het*" -> GG",
            "RSID_100" => het*" -> GG"),
        "QUERY_2" => Dict(
            "RSID_10"  => het*" -> GG",
            "RSID_100" => "AA -> GG"
        )
    )
    open(path, "w") do io
        TOML.print(io, queries)
    end
end