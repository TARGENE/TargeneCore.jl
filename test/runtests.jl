

include("ukbb.jl")
include("stackbuilding.jl")


# @testset "Test filtering_idx" begin
#     snparray = [3 2 3 0;
#                 0 2 0 1;
#                 1 1 1 1;
#                 3 0 0 1;
#                 2 3 3 0;
#                 0 0 2 2]
#     snpqueries = CSV.File(snpfile) |> DataFrame
#     snpquery = first(snpqueries)
#     snp_idx = Dict("rs3683945" => 2, "rs3706761"=> 1)
#     indices, columns = GenesInteraction.filtering_idx(snparray, snpquery, snp_idx)
#     @test columns == [2, 1]
#     @test indices == [1, 2, 4, 6]
# end

# @testset "Epistasis estimation run" begin

#     tmle_config = joinpath("config", "continuous.toml")
#     outfile = "tmleepistasis_results.csv"
#     tmleepistasis(genotypefile, 
#                     phenotypefile, 
#                     confoundersfile, 
#                     snpfile,
#                     tmle_config,
#                     outfile)
    
#     try
#         results = CSV.File(outfile) |> DataFrame
#         @test results[!, "Estimate"] isa Vector
#         @test results[!, "Standard Error"] isa Vector
#     finally
#         rm(outfile)
#     end

# end