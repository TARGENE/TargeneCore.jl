logit(X) = log.(X ./ (1 .- X))
expit(X) = 1 ./ (1 .+ exp.(-X))

"""
Hack into GLM to compute deviance on y a real
"""
function GLM.devresid(::Bernoulli, y::Vector{<:Real}, μ::Real)
    return -2*(y*log(μ) + (1-y)*log1p(-μ))
end

"""
Remove default check for y to be binary
"""
GLM.checky(y, d::Bernoulli) = nothing


function combinetotable(hot_mach::Machine, T, W)
    thot = MLJ.transform(hot_mach, T)
    return merge(thot, columntable(W))
end