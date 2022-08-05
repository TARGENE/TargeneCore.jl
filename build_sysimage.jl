using PackageCompiler
# Designed to be launched from the test directory
create_sysimage(["TargeneCore"]; 
                 sysimage_path="TargeneCoreSysimage.so",
                 precompile_execution_file="example_script.jl")

