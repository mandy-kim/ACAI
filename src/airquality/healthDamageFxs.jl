"""
# getVSL_yearly

### Inputs
- VSL_US_yearly: scalar value, US VSL value in tmoney for value year = 2020
- healthconsts: contains info for VSLtype and country-specific GDP PPP ratios (GDP PPP country_every year / GDP PPP US_2020)

### Output 
- VSL_yearly: array (number of countries) x mortality_years. VSL_yearly[i,j] = VSL of country i in mortality year j
"""
function getVSL_yearly(VSL_US_2020::Float64, ε::Float64, healthconsts::HealthConstants, simulation::Simulations)
    VSLtype = healthconsts.VSLtype #"country specific", "US", or "global average"
    GDP_PPP_ratios = healthconsts.GDP_PPP_ratios #n countries x n mortality years
    VSL_yearly = VSL_US_2020 .* (GDP_PPP_ratios .^ ε) 

    if VSLtype == "global average"
        VSL_yearly = sum(healthconsts.pop_totals .* VSL_yearly, dims=1) ./ sum(healthconsts.pop_totals) #global avg for every year
        VSL_yearly = repeat(VSL_yearly, length(healthconsts.pop_totals))
    end
    
    return VSL_yearly
end

"""
# healthdamages - translate yearly damages from mortalities for each country. If cessation lag is needed, input should already include this.

### Inputs
- ΔM_yearly: yearly mortalities for every country, including cessation lag if included
- VSL_1990_1990: US VSL for value year = 1990, in 1990-USD, sampled from distribution
- ε: income elasticity
- tmoney: damages in tmoney-year USD
- healthconsts

### Output
- Array of size ncountries x nyears - same size as input ΔM_yearly
"""
function healthdamages(i::Int64, ΔM_yearly::Array{Float64, 2}, distribution::MCDistributions, simulation::Simulations, healthconsts::HealthConstants)
    VSL_1990_1990 = sample(distribution.VSL_US_1990,i) * 1e-3 #billion dollars (distribution is in millions, conver to billions for consistency with climate damages)
    ε = sample(distribution.income_elasticity, i)
    tmoney = simulation.tmoney

    #GDP deflator: GDPdeflator_tmoney / GDPdeflator_1990
    GDP_deflator_ratio = healthconsts.GDP_deflator_ratio

    #US GDP pc ratio: GDP_US_2020 / GDP_US_1990
    US_GDP_pc_ratio = healthconsts.US_GDP_pc_ratio

    #VSL_US for 2020 in tmoney dollars -- SCALAR value
    VSL_US_2020 = VSL_1990_1990 .* (US_GDP_pc_ratio .^ ε) .* GDP_deflator_ratio
    
    #get VSL_yearly depending on VSL type
    VSL_yearly = getVSL_yearly(VSL_US_2020, ε, healthconsts, simulation)

    #calculate damages = ΔM * VSL
    damages_yearly = ΔM_yearly .* VSL_yearly #in tmoney-year USD dollars
    
    return damages_yearly
end

function discountedCost_ΔM(damages_yearly::Array{Float64,2},simulation::Simulations)
    DR = simulation.DR
    nyears = size(damages_yearly)[2]
    discount_repeat = DR .* ones(Float64, (1,nyears-1))
    discount_vector = hcat([1], exp.(-1 .* cumsum(discount_repeat, dims=2)))
    
    damages_yearly_discounted = discount_vector .* damages_yearly
    
    return damages_yearly_discounted
end