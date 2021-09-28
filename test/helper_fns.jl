
genotypefile = joinpath("data", "mouse")
phenotypefile = joinpath("data", "pheno_10_causals.txt")
confoundersfile = joinpath("data", "confouders.csv")
queryfile = joinpath("config", "query.toml")


function build_query_file(;path=queryfile, threshold=0.9)
    queries = Dict(
        "threshold" => threshold,
        "SNPS" => Dict(
            "RSID_10"  => BGEN.datadir("example.8bits.bgen"),
            "RSID_100" => BGEN.datadir("example.8bits.bgen")),
        "QUERY_1" => Dict(
            "RSID_10"  => "AG -> GG",
            "RSID_100" => "AA -> GG"),
        "QUERY_2" => Dict(
            "RSID_10"  => "AA -> GG",
            "RSID_100" => "AA -> GG"
        )
    )
    open(path, "w") do io
        TOML.print(io, queries)
    end
end