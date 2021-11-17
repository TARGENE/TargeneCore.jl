FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=Europe/Amsterdam

# Install Julia

RUN apt-get update && apt-get install -y wget

RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.3-linux-x86_64.tar.gz

RUN tar zxvf julia-1.6.3-linux-x86_64.tar.gz

RUN ln -s /julia-1.6.3/bin/julia /usr/local/bin/julia

# Installing R and hal9001
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

# Installing python3 and scikit-learn
RUN apt-get install -y python3-pip

RUN pip install scikit-learn

#Import project, build and precompile
COPY . /GenesInteraction.jl 

WORKDIR /GenesInteraction.jl

RUN julia -q --project -e'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Build the sysimage for faster startup
WORKDIR /GenesInteraction.jl/test

RUN julia -q --project=.. ../build_sysimage.jl

WORKDIR /GenesInteraction.jl

RUN alias julia='julia -Jtest/GenesInteractionSysimage.so'

RUN julia -q --project -e"using GenesInteraction"