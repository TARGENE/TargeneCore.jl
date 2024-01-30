FROM julia:1.9.0-bullseye

ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=Europe/Amsterdam

ENV JULIA_DEPOT_PATH=/opt

# Install Compiler and other tools

RUN apt-get update && \
    apt-get -y install gcc mono-mcs vim procps wget unzip lsb-release software-properties-common gnupg && \
    rm -rf /var/lib/apt/lists/*

RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"

# Install GCTA

RUN wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip

RUN unzip gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip

RUN mv gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static /usr/local/bin/gcta64

RUN chmod a+x /usr/local/bin/gcta64

# Import project and build

COPY . /TargeneCore.jl 

WORKDIR /TargeneCore.jl

RUN julia -q --project -e'using Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile(); using TargeneCore'

RUN julia --project --startup-file=no deps/build_sysimage.jl

# Check scripts are callable

RUN julia --project --startup-file=no -JTargeneCoreSysimage.so bin/prepare_confounders.jl --help

RUN julia --project --startup-file=no -JTargeneCoreSysimage.so bin/generate_tl_inputs.jl --help

RUN julia --project --startup-file=no -JTargeneCoreSysimage.so bin/generate_random_variant_estimands.jl --help
