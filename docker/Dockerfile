FROM julia:1.6.3

COPY . /GenesInteraction.jl 

WORKDIR /GenesInteraction.jl 

RUN julia --project -e"using Pkg; Pkg.instantiate(); Pkg.test()"