FROM julia:1.9.0-bullseye

ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=Europe/Amsterdam

ENV JULIA_DEPOT_PATH=/opt

RUN apt-get update && apt-get install -y wget unzip procps

# Install gcta

RUN wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip

RUN unzip gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip

RUN mv gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static /usr/local/bin/gcta64

RUN chmod a+x /usr/local/bin/gcta64

# Import project, build and precompile
COPY . /TargeneCore.jl 

WORKDIR /TargeneCore.jl

RUN julia -q --project -e'using Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile(); using TargeneCore'

RUN julia --project --startup-file=no deps/build_sysimage.jl

# Check scripts are callable 
RUN julia --project --startup-file=no -JTargeneCoreSysimage.so bin/prepare_confounders.jl --help
RUN julia --project --startup-file=no -JTargeneCoreSysimage.so bin/generate_tl_inputs.jl --help
RUN julia --project --startup-file=no -JTargeneCoreSysimage.so bin/generate_random_variant_estimands.jl --help
