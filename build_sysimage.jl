using PackageCompiler
# Designed to be launched from the test directory
create_sysimage(["TMLEEpistasis"]; 
                 sysimage_path="TMLEEpistasisSysimage.so",
                 precompile_execution_file="example_script.jl",
                 cpu_target=PackageCompiler.default_app_cpu_target())

