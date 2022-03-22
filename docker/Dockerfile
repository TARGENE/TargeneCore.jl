FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=Europe/Amsterdam

ENV JULIA_DEPOT_PATH=/opt

RUN apt-get update && apt-get install -y wget unzip

# Install gcta

RUN wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip

RUN unzip gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip

RUN mv gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static /usr/local/bin/gcta64

RUN chmod a+x /usr/local/bin/gcta64

# Install Julia

RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.2-linux-x86_64.tar.gz

RUN tar zxvf julia-1.7.2-linux-x86_64.tar.gz

RUN ln -s /julia-1.7.2/bin/julia /usr/local/bin/julia

# Install R and hal9001
RUN apt-get install -y r-base \
                    r-base-core \
                    r-recommended \
                    r-base-dev

RUN apt-get install -y libssl-dev \
                libxml2-dev \
                libcurl4-openssl-dev \
                libgit2-dev 

RUN R -e "install.packages('devtools', repos='http://cran.us.r-project.org'); \
         require(devtools);\
         install_version('hal9001', version = '0.4.1', repos = 'http://cran.us.r-project.org')"

#Import project, build and precompile
COPY . /TMLEEpistasis.jl 

WORKDIR /TMLEEpistasis.jl

RUN julia -q --project -e'using Pkg; Pkg.instantiate(); Pkg.precompile()'
