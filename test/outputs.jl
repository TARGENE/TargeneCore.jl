module TestMakeOutputs

using TargeneCore
using Test
using CairoMakie
using JLD2
using YAML

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "testutils.jl"))

@testset "Test make_outputs" begin
    tmpdir = mktempdir()
    input_prefix = joinpath(tmpdir, "tmle_output")
    output_prefix = tmpdir
    make_fake_outputs(; prefix=input_prefix)
    copy!(ARGS, [
        "make-outputs",
        input_prefix,
        string("--output-prefix=", output_prefix),
        "--verbosity=0"
    ])
    TargeneCore.julia_main()
    # Check QQ exists
    @test isfile(joinpath(tmpdir, "QQ.png"))
    # Check full results file
    results = jldopen(io ->io["results"], joinpath(tmpdir, "results.hdf5"))
    @test names(results) == ["TMLE", "OSE", "TMLE_PVALUE", "OSE_PVALUE"]
    @test size(results, 1) == 5
    # Check summary file
    summary = YAML.load_file(joinpath(tmpdir, "results.summary.yaml"))
    ## Simple Effect
    simple_effect = summary[1]
    @test simple_effect["EFFECT_TYPE"] == "AIE"
    @test simple_effect["OUTCOME"] == "High light scatter reticulocyte percentage"
    @test simple_effect["TREATMENTS"] == Dict(
            "rs10043934" => "GG => GA", 
            "RSID_103" => "GG => GA"
    )
    ose = simple_effect["OSE"]
    @test ose["PVALUE"] ≈ 3.16919e-24 atol=1e-5
    @test ose["EFFECT_SIZE"] == -1
    tmle = simple_effect["TMLE"]
    @test tmle["PVALUE"] ≈ 3.16919e-24 atol=1e-5
    @test tmle["EFFECT_SIZE"] == -1
    ## Joint Effect
    joint_effect = summary[3]
    @test joint_effect["EFFECT_TYPE"] == "AIE"
    @test joint_effect["OUTCOME"] == "High light scatter reticulocyte percentage"
    @test joint_effect["TREATMENTS"] == [
        Dict(
            "rs10043934" => "GG => GA", 
            "RSID_103" => "GG => GA"
        ),
        Dict(
            "rs10043934" => "GG => GA", 
            "RSID_103" => "GA => AA"
        )
    ]
    # These results are made up so TMLE and OSE are equal here
    for estimator in ("OSE", "TMLE")
        ose = joint_effect[estimator]
        @test ose["PVALUE"] ≈ 5.25721e-11 atol=1e-5
        effect_change_1 = ose["COMPONENTS"][1]
        @test effect_change_1["EFFECT_SIZE"] == -1.0
        @test effect_change_1["PVALUE"] ≈ 3.16919e-24 atol=1e-5
        effect_change_2 = ose["COMPONENTS"][2]
        @test effect_change_2["EFFECT_SIZE"] == -0.003
        @test effect_change_2["PVALUE"] ≈ 0.011508 atol=1e-5
    end
    ## Failed Effect
    failed_estimate = summary[5]
    @test failed_estimate["EFFECT_TYPE"] == "ATE"
    @test failed_estimate["OUTCOME"] == "L50-L54 Urticaria and erythema"
    @test failed_estimate["TREATMENTS"] == Dict(
        "RSID_104" => "GG => GA",
        "rs117913124" => "GG => GA"
    )
    ose = failed_estimate["OSE"]
    @test ose["PVALUE"] === NaN
    @test ose["EFFECT_SIZE"] == "Failed"
    tmle = failed_estimate["TMLE"]
    @test tmle["PVALUE"] === NaN
    @test tmle["EFFECT_SIZE"] == "Failed"
end

end

true