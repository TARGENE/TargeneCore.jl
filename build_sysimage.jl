using PackageCompiler
# Designed to be launched from the test directory
create_sysimage(["GenesInteraction"]; 
                 sysimage_path="GenesInteractionSysimage.so",
                 precompile_execution_file="../example_script.jl")

