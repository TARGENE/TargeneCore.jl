module TestLocoGwasEstimands

using Test
using SnpArrays
using TargeneCore
using Arrow
using DataFrames
using Serialization
using TMLE

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "tl_inputs", "test_utils.jl"))

function summary_stats_df(estimands)
    estimands = DataFrame(ESTIMAND=estimands)
    estimands.ESTIMAND_TYPE = [string(typeof(first(x.args))) for x in estimands.ESTIMAND]
    estimands.ORDER = [length(first(x.args).treatment_values) for x in estimands.ESTIMAND]
    return sort(DataFrames.combine(groupby(estimands, [:ESTIMAND_TYPE, :ORDER]), nrow), [:ORDER, :ESTIMAND_TYPE])
end
@testset "Test loco-gwas from flat list: no positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 5,
        "positivity-constraint" => 0.0,

        "%COMMAND%" => "loco-gwas", 

        "loco-gwas" => Dict{String, Any}(
            "config" => joinpath(TESTDIR, "data", "config_gwas.yaml"), 
            "traits" => joinpath(TESTDIR, "data", "traits_1.csv"),
            "pcs" => joinpath(TESTDIR, "data", "pcs.csv"),
            "bed-file" => joinpath(TESTDIR, "data", "ukbb", "loco_gwas_genotypes" ,"test_chr12.bed"), 
            "bim-file" => joinpath(TESTDIR, "data", "ukbb", "loco_gwas_genotypes" ,"test_chr12.bim"),
            "fam-file" => joinpath(TESTDIR, "data", "ukbb", "loco_gwas_genotypes" ,"test_chr12.fam")
        ),
    )
    tl_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test sort(names(trait_data)) == sort([
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", 
        "COV_1", "21003", "22001", "TREAT_1", "PC1", "PC2",
        "RSID_2","RSID_3","RSID_4","RSID_5","RSID_6","RSID_7","RSID_8","RSID_9","RSID_10","RSID_11","RSID_12","RSID_13","RSID_14","RSID_15","RSID_16","RSID_17","RSID_18","RSID_19","RSID_20","RSID_21","RSID_22","RSID_23","RSID_24","RSID_25","RSID_26","RSID_27","RSID_28","RSID_29","RSID_30","RSID_31","RSID_32","RSID_33","RSID_34","RSID_35","RSID_36","RSID_37","RSID_38","RSID_39","RSID_40","RSID_41","RSID_42","RSID_43","RSID_44","RSID_45","RSID_46","RSID_47","RSID_48","RSID_49","RSID_50","RSID_51","RSID_52","RSID_53","RSID_54","RSID_55","RSID_56","RSID_57","RSID_58","RSID_59","RSID_60","RSID_61","RSID_62","RSID_63","RSID_64","RSID_65","RSID_66","RSID_67","RSID_68","RSID_69","RSID_70","RSID_71","RSID_72","RSID_73","RSID_74","RSID_75","RSID_76","RSID_77","RSID_78","RSID_79","RSID_80","RSID_81","RSID_82","RSID_83","RSID_84","RSID_85","RSID_86","RSID_87","RSID_88","RSID_89","RSID_90","RSID_91","RSID_92","RSID_93","RSID_94","RSID_95","RSID_96","RSID_97","RSID_98","RSID_99","RSID_100","RSID_101","RSID_102","RSID_103","RSID_104","RSID_105","RSID_106","RSID_107","RSID_108","RSID_109","RSID_110","RSID_111","RSID_112","RSID_113","RSID_114","RSID_115","RSID_116","RSID_117","RSID_118","RSID_119","RSID_120","RSID_121","RSID_122","RSID_123","RSID_124","RSID_125","RSID_126","RSID_127","RSID_128","RSID_129","RSID_130","RSID_131","RSID_132","RSID_133","RSID_134","RSID_135","RSID_136","RSID_137","RSID_138","RSID_139","RSID_140","RSID_141","RSID_142","RSID_143","RSID_144","RSID_145","RSID_146","RSID_147","RSID_148","RSID_149","RSID_150","RSID_151","RSID_152","RSID_153","RSID_154","RSID_155","RSID_156","RSID_157","RSID_158","RSID_159","RSID_160","RSID_161","RSID_162","RSID_163","RSID_164","RSID_165","RSID_166","RSID_167","RSID_168","RSID_169","RSID_170","RSID_171","RSID_172","RSID_173","RSID_174","RSID_175","RSID_176","RSID_177","RSID_178","RSID_179","RSID_180","RSID_181","RSID_182","RSID_183","RSID_184","RSID_185","RSID_186","RSID_187","RSID_188","RSID_189","RSID_190","RSID_191","RSID_192","RSID_193","RSID_194","RSID_195","RSID_196","RSID_197","RSID_198","RSID_199","RSID_200"])
    @test size(trait_data) == (490, 210)
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa ComposedEstimand for e in estimands)
    summary_stats = summary_stats_df(estimands)
    @test summary_stats == DataFrame(
        ESTIMAND_TYPE=["TMLE.StatisticalATE"],
        ORDER=[1],
        nrow=[776]
    )
end

@testset "Test loco-gwas from flat list: positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 5,
        "positivity-constraint" => 0.01,

        "%COMMAND%" => "loco-gwas", 

        "loco-gwas" => Dict{String, Any}(
            "config" => joinpath(TESTDIR, "data", "config_gwas.yaml"), 
            "traits" => joinpath(TESTDIR, "data", "traits_1.csv"),
            "pcs" => joinpath(TESTDIR, "data", "pcs.csv"),
            "bed-file" => joinpath(TESTDIR, "data", "ukbb", "loco_gwas_genotypes" ,"test_chr12.bed"), 
            "bim-file" => joinpath(TESTDIR, "data", "ukbb", "loco_gwas_genotypes" ,"test_chr12.bim"),
            "fam-file" => joinpath(TESTDIR, "data", "ukbb", "loco_gwas_genotypes" ,"test_chr12.fam")
        ),
    )
    tl_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test sort(names(trait_data)) == sort([
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", 
        "COV_1", "21003", "22001", "TREAT_1", "PC1", "PC2",
        "RSID_2","RSID_3","RSID_4","RSID_5","RSID_6","RSID_7","RSID_8","RSID_9","RSID_10","RSID_11","RSID_12","RSID_13","RSID_14","RSID_15","RSID_16","RSID_17","RSID_18","RSID_19","RSID_20","RSID_21","RSID_22","RSID_23","RSID_24","RSID_25","RSID_26","RSID_27","RSID_28","RSID_29","RSID_30","RSID_31","RSID_32","RSID_33","RSID_34","RSID_35","RSID_36","RSID_37","RSID_38","RSID_39","RSID_40","RSID_41","RSID_42","RSID_43","RSID_44","RSID_45","RSID_46","RSID_47","RSID_48","RSID_49","RSID_50","RSID_51","RSID_52","RSID_53","RSID_54","RSID_55","RSID_56","RSID_57","RSID_58","RSID_59","RSID_60","RSID_61","RSID_62","RSID_63","RSID_64","RSID_65","RSID_66","RSID_67","RSID_68","RSID_69","RSID_70","RSID_71","RSID_72","RSID_73","RSID_74","RSID_75","RSID_76","RSID_77","RSID_78","RSID_79","RSID_80","RSID_81","RSID_82","RSID_83","RSID_84","RSID_85","RSID_86","RSID_87","RSID_88","RSID_89","RSID_90","RSID_91","RSID_92","RSID_93","RSID_94","RSID_95","RSID_96","RSID_97","RSID_98","RSID_99","RSID_100","RSID_101","RSID_102","RSID_103","RSID_104","RSID_105","RSID_106","RSID_107","RSID_108","RSID_109","RSID_110","RSID_111","RSID_112","RSID_113","RSID_114","RSID_115","RSID_116","RSID_117","RSID_118","RSID_119","RSID_120","RSID_121","RSID_122","RSID_123","RSID_124","RSID_125","RSID_126","RSID_127","RSID_128","RSID_129","RSID_130","RSID_131","RSID_132","RSID_133","RSID_134","RSID_135","RSID_136","RSID_137","RSID_138","RSID_139","RSID_140","RSID_141","RSID_142","RSID_143","RSID_144","RSID_145","RSID_146","RSID_147","RSID_148","RSID_149","RSID_150","RSID_151","RSID_152","RSID_153","RSID_154","RSID_155","RSID_156","RSID_157","RSID_158","RSID_159","RSID_160","RSID_161","RSID_162","RSID_163","RSID_164","RSID_165","RSID_166","RSID_167","RSID_168","RSID_169","RSID_170","RSID_171","RSID_172","RSID_173","RSID_174","RSID_175","RSID_176","RSID_177","RSID_178","RSID_179","RSID_180","RSID_181","RSID_182","RSID_183","RSID_184","RSID_185","RSID_186","RSID_187","RSID_188","RSID_189","RSID_190","RSID_191","RSID_192","RSID_193","RSID_194","RSID_195","RSID_196","RSID_197","RSID_198","RSID_199","RSID_200"])
    @test size(trait_data) == (490, 210)
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa ComposedEstimand for e in estimands)
    summary_stats = summary_stats_df(estimands)
    @test summary_stats == DataFrame(
        ESTIMAND_TYPE=["TMLE.StatisticalATE"],
        ORDER=[1],
        nrow=[616]
    )
end

end

true