FROM julia:1.8.2-bullseye

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

RUN julia -q --project=/TargeneCore.jl -e'using Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile(); using TargeneCore'

RUN julia -q --project=/TargeneCore.jl -e'using Pkg; Pkg.test()'

