abstract type RandomVariable end

## Discrete Random Variables
#Constant Distribution
struct ConstantDist <: RandomVariable
    value::Float64
end

function get_distribution(x::ConstantDist)
    return DiscreteNonParametric([x.value],[1])
end

#Discrete Distribution
struct DiscreteDist <: RandomVariable
    S::Array{Float64}
    a::Int64
    b::Int64
end

function get_distribution(x::DiscreteDist)
    return DiscreteUniform(x.a,x.b)
end


## Continuous Random Variables
#Uniform Distribution
struct UniformDist <: RandomVariable
    S::Array{Float64}
    a::Float64
    b::Float64
end

function get_distribution(x::UniformDist)
    return Uniform(x.a,x.b)
end

#Triangular Distribution
struct TriangDist <: RandomVariable
    S::Array{Float64}
    min::Float64
    top::Float64
    max::Float64
end

function get_distribution(x::TriangDist)
    return TriangularDist(x.min, x.max, x.top)
end

#Normal Distribution
struct NormalDist <: RandomVariable
    S::Array{Float64}
    mu::Float64
    sigma::Float64
    lower_bound::Float64
    upper_bound::Float64
end

function get_distribution(x::NormalDist)
    temp_dist = Normal(x.mu, x.sigma)
    return truncated(temp_dist, x.lower_bound, x.upper_bound)
end

#Weibull Distribution
struct WeibullDist <: RandomVariable
    S::Array{Float64}
    shape::Float64
    scale::Float64
end

function get_distribution(x::WeibullDist)
    return Weibull(x.shape,x.scale)
end

#RoeBaker Distribution
struct RoeBakerDist <: RandomVariable
    S::Array{Float64}
    f_mu::Float64
    f_sigma::Float64
    RoeB_max::Float64
end


#"random" sampling functions
function sample(x::ConstantDist,i)
    return x.value
end

function sample(x::RandomVariable, i)
	return quantile(get_distribution(x), x.S[i])
end

function sample(x::RoeBakerDist, i)
    f_dist = Normal(x.f_mu, x.f_sigma)
    f_max = 1 - 1.2/x.RoeB_max
    f_dist = truncated(f_dist, -Inf, f_max)

    return 1.2/(1-quantile(f_dist,x.S[i]))
end


function generate_RVs(dict::Dict, S::AbstractArray)
    distribution_type = dict["distribution"]
    if distribution_type == "Constant"
        a = dict["value"]
        return ConstantDist(a)
        
    elseif distribution_type == "Uniform"
        a = dict["a"]
        b = dict["b"]
        return UniformDist(S, a, b)
        
    elseif distribution_type == "Triangular"
        min = dict["min"]
        top = dict["top"]
        max = dict["max"]
        return TriangDist(S, min, top, max)
        
    elseif distribution_type == "DiscreteUniform"
        a = dict["a"]
        b = dict["b"]
        return DiscreteDist(S, a, b)
        
    elseif distribution_type == "RoeBaker"
        f_mu = dict["f_mu"]
        f_sigma = dict["f_sigma"]
        max_val = dict["RoeB_max"]
        return RoeBakerDist(S, f_mu, f_sigma, max_val)

    elseif distribution_type == "Weibull"
        shape = dict["shape"]
        scale = dict["scale"]
        return WeibullDist(S, shape,scale)
        
    elseif distribution_type == "Normal"
        mu = dict["mu"]
        sigma = dict["sigma"]
        lb = dict["lower_bound"]
        ub = dict["upper_bound"]
        if lb == "-Inf"
            lb = -Inf
        end
        if ub == "Inf"
            ub = Inf
        end
        return NormalDist(S, mu, sigma, lb, ub)
    else 
        error("Unknown distribution")
    end
end
