FROM julia:1.6.3

COPY . /GenesInteraction.jl 

WORKDIR /GenesInteraction.jl

RUN julia -q --project -e'using Pkg; Pkg.instantiate(); Pkg.precompile()'

WORKDIR /GenesInteraction.jl/test

RUN julia -q --project=.. ../build_sysimage.jl

WORKDIR /GenesInteraction.jl

RUN alias julia='julia -Jtest/GenesInteractionSysimage.so'

RUN julia -e"using GenesInteraction"