module TestGRM
using Test
using TargeneCore

@testset "Test readGRM" begin
    prefix = joinpath("data", "grm", "test.grm")
    GRM, ids = TargeneCore.readGRM(prefix)
    @test size(GRM, 1) == 18915
    @test size(ids, 1) == 194
end

end;

true