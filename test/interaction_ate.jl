
using Test
using MLJ


@testset "Test tomultivariate" begin
    T = (t1 = categorical([true, false, false, true]),
         t2 = categorical([true, true, false, false]))
    # The mapping is fixed and harcoded
    @test GenesInteraction.tomultivariate(T) == categorical([1, 3, 4, 2])
end

@testset "Testing the various intermediate functions" begin
    T = (t1 = categorical([true, false, false, true, true, true, false]),
         t2 = categorical([true, true, false, false, false, true, false]))
    y = categorical([true, true, false, false, false, true, true])
    W = (W1=categorical([true, true, false, true, true, true, true]),)
    
    dummy_model = @pipeline OneHotEncoder() ConstantClassifier()
end