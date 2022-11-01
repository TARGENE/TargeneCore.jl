FROM julia:1.8.2-bullseye

ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=Europe/Amsterdam

ENV JULIA_DEPOT_PATH=/opt

RUN apt-get update && apt-get install -y wget unzip procps libcurl4-openssl-dev libbz2-dev liblzma-dev zlib1g

# Install gcta

RUN wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip

RUN unzip gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip

RUN mv gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static /usr/local/bin/gcta64

RUN chmod a+x /usr/local/bin/gcta64

# Install BGEN


RUN wget http://code.enkre.net/bgen/tarball/release/bgen.tgz

RUN tar -xzf bgen.tgz

RUN cd bgen.tgz && \
    ./waf configure && \
    ./waf && \
    ./build/test/test_bgen && \
    ./waf install

# Install R and hal9001

RUN apt-get install -y r-base \
                    r-base-core \
                    r-recommended \
                    r-base-dev

RUN apt-get install -y libssl-dev \
                libxml2-dev \
                libcurl4-openssl-dev \
                libgit2-dev \
                libfontconfig1-dev \
                libharfbuzz-dev \
                libfribidi-dev \
                libpng-dev \
                libtiff5-dev \
                libjpeg-dev

RUN R -e "install.packages('devtools', repos='http://cran.us.r-project.org'); \
         require(devtools);\
         install_version('hal9001', version = '0.4.1', repos = 'http://cran.us.r-project.org')"

# Import project, build and precompile

COPY . /TargeneCore.jl 

RUN julia -q --project=/TargeneCore.jl -e'using Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile(); using TargeneCore'

# Prepend a writable DEPOT_PATH ($PWD is mounted by default by singularity)

RUN mkdir $PWD/julia_depot

ENV JULIA_DEPOT_PATH=$PWD/julia_depot:$JULIA_DEPOT_PATH

# RUN julia -q --project build_sysimage.jl

# RUN echo 'alias julia="julia -J TMLEEpistasisSysimage.so"' >> ~/.bashrc
