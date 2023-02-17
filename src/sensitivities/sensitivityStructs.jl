Base.@kwdef struct TotalImpacts
    multipliers::Array{Float64}
    RF::Array{Float64,1}
    PM25::Array{Float64,3}
    Sfc_O3::Array{Float64,3}
    Col_O3::Array{Float64,3}
    Zonal_O3::Array{Float64,3}
end

Base.@kwdef struct SensitivityRegionImpacts 
    RF::Array{Float64,1}
    PM25::Array{Float64,3}
    Sfc_O3::Array{Float64,3}
    Col_O3::Array{Float64,3}
    Zonal_O3::Array{Float64,3}
    mask::Array{Float64,2}
    lev_47::Array{Int64,1}
    lev_72::Array{Int64,1}
end

Base.@kwdef struct EmissionSensitivities
    n_case::Int64
    sensitivity_dict::Dict{Int64,SensitivityRegionImpacts}
end


# IF sensitivities give step-response impact results - the following functions calculates the impulse response as the Δimpact between each year. 
function calculateImpulse(RF::Array{Float64,1})
    n_year = length(RF)
    impact_impulse = cat(RF[1],cat([RF[i].-RF[i-1] for i in 2:n_year]...;dims=1);dims=1)
    return impact_impulse
end

function calculateImpulse(impact::Array{Float64,3})
    n_year = size(impact,3)
    impact_impulse = cat(impact[:,:,1],cat([impact[:,:,i].-impact[:,:,i-1] for i in 2:n_year]...;dims=3);dims=3)
    return impact_impulse
end