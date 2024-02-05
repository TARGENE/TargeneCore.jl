using PackageCompiler
PackageCompiler.create_sysimage(
    ["TargeneCore"], 
    cpu_target="generic",
    sysimage_path="TargeneCoreSysimage.so", 
    precompile_execution_file="deps/execute.jl", 
)
