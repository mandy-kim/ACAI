Base.@kwdef struct Emissions_orig
    FB::Array{Float64,1} #kg
    CO2::Array{Float64,1} #kg
    NOx::Array{Float64,1} #kg
    lc_CO2::Array{Float64,1} #kg
    lc_N2O::Array{Float64,1} #kg
    lc_CH4::Array{Float64,1} #kg
    lc_CH4_ff::Array{Float64,1} #kg
end

function getYearlyEmis(simulation::Simulations_orig)
    years = simulation.emissionyears
    totalimpact = simulation.totalimpact
    FB_array = zeros(totalimpact)
    CO2_array = zeros(totalimpact)
    NOx_array = zeros(totalimpact)
    lc_CO2_array = zeros(totalimpact)
    lc_N2O_array = zeros(totalimpact)
    lc_CH4_array = zeros(totalimpact)
    lc_CH4_ff_array = zeros(totalimpact)

    ni = 0
    for yeari in years
        ni+=1
        FB_array[ni] = simulation.FB_inv[ni]
        CO2_array[ni] = simulation.CO2_inv[ni]
        NOx_array[ni] = simulation.NOx_inv[ni]
        lc_CO2_array[ni] = simulation.lc_CO2_inv[ni]
        lc_N2O_array[ni] = simulation.lc_N2O_inv[ni]
        lc_CH4_array[ni] = simulation.lc_CH4_inv[ni]
        lc_CH4_ff_array[ni] = simulation.lc_CH4_inv[ni] .* simulation.lc_frac_CH4_ff[ni]
    end

    emission = Emissions_orig(
        FB = FB_array, 
        CO2 = CO2_array, 
        NOx = NOx_array,
        lc_CO2 = lc_CO2_array,
        lc_N2O = lc_N2O_array,
        lc_CH4 = lc_CH4_array,
        lc_CH4_ff = lc_CH4_ff_array,
        )
    return emission
end