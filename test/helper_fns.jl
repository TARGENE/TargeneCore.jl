
genotypefile = joinpath("data", "mouse")
phenotypefile = joinpath("data", "pheno_10_causals.txt")
confoundersfile = joinpath("data", "confouders.csv")
queryfile = joinpath("data", "query.csv")


function build_query_file(path)
    CSV.write(path, DataFrame(
            RSID     = ["RSID_10", "RSID_100"],
            BGENFILE = [BGEN.datadir("example.8bits.bgen"), BGEN.datadir("example.8bits.bgen")],
            ALLELE_1 = ["A/G", "G/G"],
            ALLELE_2 = ["A/A", "A/A"]
            ))
end